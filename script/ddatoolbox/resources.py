#-------------------------------------------------------------------------------
# Name:        resources.py
# Purpose:     Loaders for cached reference resources used across the pipeline:
#              PubChem compound-property lookups (live and cached), non-natural-
#              product and InChIKey-to-ClassyFire lookup tables, and the bundled
#              compound name list.
#
# Author:      Fumio Matsuda@ShinMassBank PJ
#
# Created:     23/07/2023
# Copyright:   (c) hinMassBank PJ 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from pathlib import Path
import glob, csv, pickle
import pubchempy as pcp
from . import ddatoolbox
from . import paths
from importlib.resources import read_text


def get_pubchemdatafrominchikey(inchikey):
    """Fetch compound properties from PubChem for a given InChIKey (live lookup).

    Queries PubChem (via ``pubchempy``) for compounds matching ``inchikey``, and if
    at least one match is found, fetches the full record of the lowest-CID compound
    and extracts a normalized subset of properties. Falls back to scanning the raw
    PubChem record's ``props`` list for ``SMILES``/``IUPAC Name`` when the
    convenience attributes are absent from the compound dict (older PubChemPy
    versions).

    Args:
        inchikey (str): InChIKey of the compound to look up.

    Returns:
        dict: A dictionary with keys ``"ExactMass"``, ``"MolecularFormula"``,
            ``"InChI"``, ``"InChIKey"``, ``"CID"``, ``"IsomericSMILES"``, and
            ``"IUPACName"`` (the latter two only when resolvable). Returns ``{}``
            if no matching compound is found on PubChem.
    """
    pdcompoundall = pcp.get_compounds(inchikey, 'inchikey')
    cids = sorted([int(x.cid) for x in pdcompoundall])
    if len(cids) >= 1:
        c = pcp.Compound.from_cid(cids[0])
        pdproperty_temp = c.to_dict()
        pdproperty = {}
        pdproperty['ExactMass']= pdproperty_temp['exact_mass']
        pdproperty['MolecularFormula']= pdproperty_temp['molecular_formula']
        pdproperty['InChI']= pdproperty_temp['inchi']
        pdproperty['InChIKey']= pdproperty_temp['inchikey']
        pdproperty['CID']= pdproperty_temp['cid']
        if 'smiles' in pdproperty_temp:
            pdproperty['IsomericSMILES'] = pdproperty_temp['smiles']
        else:
            for item in pdproperty_temp['record']['props']:
                if item['urn']['label'] == 'SMILES' and item['urn']['name'] == 'Absolute':
                    pdproperty['IsomericSMILES'] = item['value']['sval']
        if 'iupac_name' in pdproperty_temp:
            pdproperty['IUPACName'] = pdproperty_temp['iupac_name']
        else:
            for item in pdproperty_temp['record']['props']:
                if item['urn']['label'] == 'IUPAC Name':
                    pdproperty['IUPACName'] = item['value']['sval']


    else:
        pdproperty = {}
    return(pdproperty)





def load_pubchemchecked():
    """Load the on-disk cache of previously fetched PubChem compound data.

    Reads the pickled dictionary at ``paths.generate_paths_resource()["pubchemchecked"]``,
    avoiding repeated live PubChem queries (see :func:`get_pubchemdatafrominchikey`)
    for compounds already looked up in a previous run.

    Returns:
        dict: Mapping of InChIKey to the PubChem property dictionary previously
            fetched for it (same shape as :func:`get_pubchemdatafrominchikey`'s
            return value).
    """
    respath = paths.generate_paths_resource()
    pubchem_checked = {}
    with open(respath['pubchemchecked'], 'rb') as f:
        pubchem_checked = pickle.load(f)
    return(pubchem_checked)

def load_inchikey2classyfire():
    """Load the InChIKey-to-ClassyFire direct-parent-class lookup table.

    Reads the ``inchikey2classyfire`` resource table (``InchikeyClassyfireDB``,
    via ``ddatoolbox.ddatoolbox.readsmalltable``) and builds a dictionary keyed by
    the first (14-character) block of each InChIKey, used by the pipeline's
    ``target == "class"`` annotation logic (see ``DDAFdr.py`` /
    ``DDARetrieve.py``) to map a compound hit to its ClassyFire compound class.

    Returns:
        dict: Mapping of short InChIKey (first block, str) to ClassyFire
            "Direct parent name" (str).
    """
    #
    # reading InchikeyClassyfireDB-VS5
    #
    #softwarepath = "../software/"
    #softwareroot = Path(softwarepath).resolve()
    #filenametemp = softwareroot/"InchikeyClassyfireDB-VS5.txt"

    respath = paths.generate_paths_resource()
    filenametemp = respath["inchikey2classyfire"]
    csv_db = ddatoolbox.readsmalltable(filenametemp, keywords = ["InChIKey", "Direct parent name", "Direct parent ID"], filetype = "tab")
    InChIkeydb = {}

    for tempdata in csv_db:
        InChIkey = tempdata['InChIKey']
        InChIkey_first = InChIkey.split("-")[0]
        InChIkeydb[InChIkey_first] = tempdata["Direct parent name"]

    return(InChIkeydb)
def load_compoundlist():
    """Load the bundled reference compound-name list shipped with the package.

    Reads ``compoundlist.txt`` from the ``ddatoolbox`` package data (via
    ``importlib.resources.read_text``), split into lines.

    Returns:
        list[str]: Lines of ``compoundlist.txt``, in file order.
    """
    reader = read_text("ddatoolbox", "compoundlist.txt", encoding="utf-8").split('\n')
    sensitive = []
    unrealistic = {}
    for line in reader:
        temp = line.split("\t")
        inchikey = temp[0]
        reason = temp[1]
        if reason == "sensitive":
            sensitive.append(inchikey)
        if reason == "unrealistic":
            text = temp[2]
            unrealistic[inchikey] = text 

    return(sensitive, unrealistic)

def save_pubchemchecked(pubchem_checked):
    """Persist the PubChem compound-data cache to disk as a pickle file.

    Writes ``pubchem_checked`` to
    ``paths.generate_paths_resource()["pubchemchecked"]``, overwriting any
    existing cache file. Intended to be called after merging newly fetched
    PubChem lookups (see :func:`get_pubchemdatafrominchikey`) into the dictionary
    returned by :func:`load_pubchemchecked`.

    Args:
        pubchem_checked (dict): Mapping of InChIKey to PubChem property dictionary
            to be cached.

    Returns:
        None
    """
    respath = paths.generate_paths_resource()
    with open(respath['pubchemchecked'],"wb") as f:
        pickle.dump(pubchem_checked, f)