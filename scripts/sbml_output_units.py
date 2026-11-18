"""Describe DNAr's imported state units without changing the SBML importer."""
import json
import sys

import libsbml

document = libsbml.readSBML(sys.argv[1])
model = document.getModel()
if model is None:
    raise ValueError("SBML document has no model")
species = {}
for item in model.getListOfSpecies():
    compartment = model.getCompartment(item.getCompartment())
    size = compartment.getSize() if compartment.isSetSize() else None
    species[item.getId()] = {
        "amount_state": item.getHasOnlySubstanceUnits(),
        "volume": size,
    }
print(json.dumps(species, allow_nan=False))
