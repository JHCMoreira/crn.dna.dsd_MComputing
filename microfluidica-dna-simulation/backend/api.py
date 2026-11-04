"""Local-first API. Run one API process: one durable queue, isolated workers."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from pathlib import Path
import importlib.util
import json
import os
import sqlite3
import subprocess
import sys
import signal
import threading
import uuid

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from .models import Project, import_project
from .engine import atomic_json, read_json
from .cad import compile_netlist

ROOT = Path(__file__).resolve().parents[1]
RUNS = Path(os.environ.get("FLUIDNA_RUNS", ROOT / ".runs")).resolve()
HOST_RUNS = Path(os.environ.get("FLUIDNA_HOST_RUNS", RUNS)).resolve()
RUNS.mkdir(parents=True, exist_ok=True)
DB = RUNS / "jobs.sqlite3"
pool = ThreadPoolExecutor(max_workers=1)
lock = threading.Lock()
processes = {}


def db():
    connection = sqlite3.connect(DB)
    connection.row_factory = sqlite3.Row
    return connection


with db() as conn:
    conn.execute("CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, status TEXT, error TEXT, created TEXT DEFAULT CURRENT_TIMESTAMP)")
    conn.execute("UPDATE jobs SET status='failed', error='Service restarted before completion' WHERE status IN ('queued','running')")


def get_job(jid):
    with db() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
    if row is None:
        raise HTTPException(404, "Run not found")
    return dict(row)


def status(jid, value, error=None):
    with db() as conn:
        conn.execute("UPDATE jobs SET status=?, error=? WHERE id=? AND status!='cancelled'", (value, error, jid))


def execute(jid):
    folder = RUNS / jid
    proc = None
    log = None
    docker = not local_engine()
    try:
        with lock:
            if get_job(jid)["status"] == "cancelled":
                return
            status(jid, "running")
            if docker:
                # Create before publishing the process handle. Cancellation can
                # then remove a known container even if it has not started yet.
                created = subprocess.run(["docker", "create", "--rm", "--init", "--name", f"fluidna-{jid}",
                    "--network", "none", "--memory", "2g", "--cpus", "2", "--read-only",
                    "--tmpfs", "/tmp:rw,nosuid,nodev,size=128m", "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
                    "-v", f"{HOST_RUNS / jid}:/run",
                    os.environ.get("FLUIDNA_ENGINE_IMAGE", "fluidna-engine:v2"), "python3", "-m", "backend.worker", "/run"],
                    capture_output=True, text=True, timeout=30)
                if created.returncode:
                    raise RuntimeError(created.stderr)
                cmd = ["docker", "start", "--attach", f"fluidna-{jid}"]
            else:
                cmd = [sys.executable, "-m", "backend.worker", str(folder)]
            log = (folder / "worker.log").open("w", encoding="utf-8")
            proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                    start_new_session=(os.name != "nt"))
            processes[jid] = (proc, docker)
        timeout = read_json(folder / "project.json")["settings"]["timeout"]
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            stop(proc, docker, jid)
            raise RuntimeError(f"Simulation exceeded {timeout}s; no result published")
        finally:
            log.close()
        if code != 0:
            detail = (folder / "worker.log").read_text(encoding="utf-8", errors="replace")[-5000:]
            raise RuntimeError(f"Worker failed ({code}). {detail}")
        result = read_json(folder / "result.json")
        if result.get("schema_version") != 2:
            raise RuntimeError("Invalid worker result")
        status(jid, "completed")
    except Exception as exc:
        status(jid, "failed", str(exc))
    finally:
        if log is not None:
            log.close()
        # Also clean up a created container if launching its attach process failed.
        if docker and proc is None:
            try:
                subprocess.run(["docker", "rm", "--force", f"fluidna-{jid}"],
                               capture_output=True, timeout=15)
            except (OSError, subprocess.TimeoutExpired):
                pass
        with lock:
            processes.pop(jid, None)


def local_engine():
    try:
        return importlib.util.find_spec("mmft.simulator") is not None
    except ModuleNotFoundError:
        return False


def stop(proc, docker, jid):
    if docker:
        subprocess.run(["docker", "rm", "--force", f"fluidna-{jid}"], capture_output=True, timeout=15)
    elif os.name == "nt":
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
    else:
        os.killpg(proc.pid, signal.SIGTERM)
    if proc.poll() is None:
        proc.kill()
    proc.wait(timeout=15)


@asynccontextmanager
async def lifespan(app):
    yield
    with lock:
        for jid, (proc, docker) in list(processes.items()):
            stop(proc, docker, jid)
    pool.shutdown(wait=False, cancel_futures=True)


app = FastAPI(title="fluiDNA v2", lifespan=lifespan)


@app.get("/api/v2/health")
def health():
    return {"status": "ready", "engine_transport": "local" if local_engine() else "docker", "schema_version": 2}


@app.post("/api/v2/projects/import")
def import_file(data: dict):
    try:
        return import_project(data).model_dump()
    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(422, str(exc)) from exc


@app.post("/api/v2/projects/validate")
def validate(project: Project):
    try:
        project.simulation_ready()
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    return {"valid": True}


@app.post("/api/v2/cad/compile")
def compile_cad(data: dict):
    """Compile a junction/netlist description into an editable v2 project."""
    try:
        result = compile_netlist(data)
        Project.model_validate(result["project"])
        return result
    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(422, str(exc)) from exc


@app.post("/api/v2/runs", status_code=202)
def submit(project: Project):
    validate(project)
    jid = uuid.uuid4().hex
    folder = RUNS / jid
    folder.mkdir()
    atomic_json(folder / "project.json", project.model_dump())
    with db() as conn:
        conn.execute("INSERT INTO jobs(id,status) VALUES (?, 'queued')", (jid,))
    pool.submit(execute, jid)
    return get_job(jid)


@app.get("/api/v2/runs/{jid}")
def job(jid: str):
    return get_job(jid)


@app.get("/api/v2/runs")
def history():
    with db() as conn:
        entries = [dict(row) for row in conn.execute("SELECT * FROM jobs ORDER BY created DESC, rowid DESC LIMIT 100")]
    for entry in entries:
        try:
            project = read_json(RUNS / entry["id"] / "project.json")
            entry.update(name=project["name"], mixing=project["settings"]["mixing"])
        except (OSError, ValueError, KeyError):
            entry.update(name=entry["id"], mixing="unknown")
    return entries


@app.post("/api/v2/runs/{jid}/cancel")
def cancel(jid: str):
    with lock:
        info = get_job(jid)
        if info["status"] not in ("queued", "running"):
            return info
        status(jid, "cancelled")
        if jid in processes:
            stop(*processes[jid], jid)
    return get_job(jid)


@app.get("/api/v2/runs/{jid}/result")
def result(jid: str):
    if get_job(jid)["status"] != "completed":
        raise HTTPException(409, "Run has no verified completed result")
    return FileResponse(RUNS / jid / "result.json", media_type="application/json")


@app.get("/api/v2/runs/{jid}/log")
def log(jid: str):
    get_job(jid)
    path = RUNS / jid / "worker.log"
    return {"log": path.read_text(encoding="utf-8", errors="replace")[-20000:] if path.exists() else "Queued"}


@app.get("/api/v2/examples")
def examples():
    return [p.name for p in sorted((ROOT / "examples").glob("*.json"))]


@app.get("/api/v2/examples/{name}")
def example(name: str):
    if name not in examples():
        raise HTTPException(404)
    return import_file(read_json(ROOT / "examples" / name))


if (ROOT / "frontend/studio/dist").exists():
    app.mount("/", StaticFiles(directory=ROOT / "frontend/studio/dist", html=True), name="studio")
