# Testes de compatibilidade SBML

Automação da SBML Test Suite para o Papel A: preparar o ambiente, executar
os casos semânticos e organizar os resultados por recurso do SBML.

## Como executar

Na raiz do repositório, com Docker disponível:

```powershell
docker build -t dnar-sbml -f .devcontainer/Dockerfile.sbml .
docker run --rm -v "${PWD}:/workspaces/dnar" dnar-sbml Rscript scripts/download_sbml_test_suite.R
docker run --rm -v "${PWD}:/workspaces/dnar" dnar-sbml Rscript scripts/run_sbml_baseline.R
```

No Bash, substituir `${PWD}` por `$PWD`. Com R e libSBML já instalados,
os dois scripts podem ser executados diretamente com `Rscript`.

A suíte fica em `external/sbml-test-suite`, fora do versionamento. O download
usa a versão 3.5.0, commit `cf38585fac5de8e0e90112febb62851ee2181816`.

Para conferir o ambiente com três casos:

```bash
Rscript scripts/run_sbml_baseline.R --cases=00001,00025,00048
```

## Método

Cada caso usa uma variante SBML, preferindo L3V2. Quando ela não existe,
é escolhida a versão mais recente disponível. Os tempos e variáveis vêm
do arquivo `settings.txt`; quantidade e concentração são ajustadas usando
o tamanho do compartimento fixo.

A simulação usa `react()` com deSolve/LSODA, `rtol=1e-10` e `atol=1e-12`.
A comparação segue o critério por ponto:

```text
abs(referência - resultado) <= absoluto + relativo * abs(referência)
```

Os limites da rodada são `1e-6` absoluto e `1e-4` relativo, conforme proposta
do grupo. A coluna `passed_suite_tolerances` também registra a comparação
com os limites originais de cada caso.

Regras, eventos, initialAssignment, definições de função e pacotes ainda
ausentes no importador são marcados como `nao suportado` antes da simulação.
Essa lista deve ser revista após as correções de B/C. A opção
`--no-preflight` desativa o filtro para investigar um caso.

## Resultados

A rodada inicial está em [reports/sbml-baseline](../../reports/sbml-baseline/RESUMO.md).
Novas execuções escrevem em `results/sbml/`:

| Arquivo | Conteúdo |
| --- | --- |
| `cases.csv` | Caso, nível/versão, validação, erro, status, tags e mensagem |
| `status.csv` | Contagem por resultado |
| `component_tags.csv`, `test_tags.csv` | Contagens por tags da suíte |
| `causes.csv` | Triagem das falhas por tags e mensagens |
| `compatibility.csv` | Resultados por nível/versão selecionado |
| `manifest.json`, `sessionInfo.txt` | Versões e parâmetros da execução |

Os resultados individuais ficam em `cases-rds/`. Para reconstruir as
tabelas sem repetir a simulação:

```bash
Rscript scripts/summarize_sbml_baseline.R results/sbml
```

`--output-dir=CAMINHO` permite guardar uma rodada nova em outra pasta.
Isso permite comparar a base com os próximos merges sem sobrescrever os
resultados anteriores. Em Windows/OneDrive, a leitura da suíte pode ser lenta;
uma cópia no disco Linux do container evita esse custo.

## Verificação

```powershell
docker run --rm -v "${PWD}:/workspaces/dnar" dnar-sbml Rscript -e "testthat::test_file('tests/testthat/test_sbml_testbench.R', stop_on_failure=TRUE)"
```

Os testes cobrem leitura dos casos, tolerâncias, tempos, conversão de unidades
e consolidação de resultados com erro ou recurso não suportado.

## Referências

- [Formato dos casos e tags](https://sbml.org/software/sbml-test-suite/docs/test-case-details/)
- [Critério de avaliação](https://sbml.org/software/sbml-test-suite/docs/evaluating-outputs/)
- [SBML Test Suite 3.5.0](https://github.com/sbmlteam/sbml-test-suite/releases/tag/3.5.0)
