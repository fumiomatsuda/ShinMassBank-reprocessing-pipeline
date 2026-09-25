#-------------------------------------------------------------------------------
# Name:        mfsearch.py
# Purpose:     A toolbox for molecular search and its related functioons
#
# Author:      Fumio Matsuda@ShinMassBank PJ
#
# Created:     23/07/2023
# Copyright:   (c) hinMassBank PJ 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import csv, math, re, os

#: Exact (monoisotopic) mass of each supported element, in Da.
elements = {
"Cl" : 34.968853,
"Br" : 78.918336,
"Si" : 27.976928,
"C" : 12,
"H" : 1.0078250319,
"N" : 14.0030740074,
"O" : 15.994914622,
"S" : 31.972071,
"P" : 30.97376,
"F" : 18.998403
}


def mfsearch(mz, threshold = 2, ihdtype = "non-integer", Climit = 99, Hlimit = 999, Nlimit = 99, Olimit = 99, Plimit = 3, Slimit = 0, polarity = "positive"):
    """Search candidate molecular formulas matching a given m/z value.

    Looks up the pre-generated exact-mass reference table for the nominal mass
    nearest ``mz`` (``./res/exactmass/<nominal_mass>.txt``, produced by
    ``DDAtools.exactmassgenerator``), then filters candidate formulas by IHD
    parity (integer vs. non-integer, to match even-/odd-electron ion rules),
    per-element atom-count limits, and the m/z matching ``threshold`` (in ppm).
    The reference formula's exact mass is polarity-corrected by the electron mass
    (``0.0005486``) before comparison.

    Args:
        mz (float): Target m/z value to search for.
        threshold (float, optional): Matching tolerance in ppm. Defaults to ``2``.
        ihdtype (str, optional): ``"integer"`` to keep only formulas with an
            integer index of hydrogen deficiency (IHD), or ``"non-integer"`` to
            keep only formulas with a half-integer IHD. Defaults to
            ``"non-integer"``.
        Climit (int, optional): Maximum carbon atom count. Defaults to ``99``.
        Hlimit (int, optional): Maximum hydrogen atom count. Defaults to ``999``.
        Nlimit (int, optional): Maximum nitrogen atom count. Defaults to ``99``.
        Olimit (int, optional): Maximum oxygen atom count. Defaults to ``99``.
        Plimit (int, optional): Maximum phosphorus atom count. Defaults to ``3``.
        Slimit (int, optional): Maximum sulfur atom count. Defaults to ``0``.
        polarity (str, optional): ``"POSITIVE"``, ``"NEGATIVE"``, or any other
            value for neutral mass (no electron-mass correction).
            Defaults to ``"positive"`` (note: comparison is case-sensitive, so this
            default does not actually match ``"POSITIVE"``/``"NEGATIVE"`` and
            behaves as the neutral-mass case unless the caller passes the exact
            uppercase string).

    Returns:
        tuple: A tuple of two elements.

            - result (dict): Mapping of molecular formula (str) to a dict with keys
              ``"delta"`` (signed m/z deviation), ``"adsdelta"`` (absolute m/z
              deviation, used for sorting), ``"mz_corrected"``, ``"ihd"``,
              ``"exactmass"``, and per-element atom counts ``"C"``, ``"H"``,
              ``"N"``, ``"O"``, ``"P"``, ``"S"``.
            - result_list (list[str]): Formula keys of ``result``, sorted by
              ascending absolute m/z deviation.

    Note:
        Returns ``({}, [])`` when the nominal-mass reference file for ``mz`` does
        not exist under ``./res/exactmass/``.
    """

    result = {}
    result_list = sorted(result.keys(), key=lambda x: result[x]['delta'])

    threshold_value = float(mz) * threshold /1000000
    mz_corrected = mz

    mztext = str(int( mz_corrected+0.2))
    filename = "./res/exactmass/"+mztext+".txt"
    if not os.path.isfile(filename):
        return(result, result_list)

    f = open(filename, 'r',  encoding='UTF-8')
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        temprow = list(row)
        mf = temprow[0]

        if polarity == "POSITIVE":
            exactmass = float(temprow[2]) - 0.0005486
        elif polarity == "NEGATIVE":
            exactmass = float(temprow[2]) + 0.0005486
        else:
            exactmass = float(temprow[2])

        ihd = float(temprow[1])
        if ihdtype == "non-integer":
            if (abs((ihd * 2) % 2) == 0.0):
                continue
        if ihdtype == "integer":
            if (abs((ihd * 2) % 2) == 1.0):
                continue

        Cnumber = int(temprow[3])
        if Cnumber > int(Climit):
            continue
        Hnumber = int(temprow[4])
        if Hnumber > int(Hlimit):
            continue
        Nnumber = int(temprow[5])
        if Nnumber > int(Nlimit):
            continue
        Onumber = int(temprow[6])
        if Onumber > int(Olimit):
            continue
        Pnumber = int(temprow[7])
        if Pnumber > int(Plimit):
            continue
        Snumber = int(temprow[8])
        if Snumber > int(Slimit):
            continue
        delta = abs(mz_corrected - exactmass)
        if (delta < threshold_value):
            result[mf] = {'delta':  mz_corrected-exactmass,'adsdelta':  delta, 'mz_corrected': mz_corrected,  'ihd':ihd, 'exactmass':exactmass, 'C':Cnumber, 'H':Hnumber, 'N':Nnumber, 'O':Onumber, 'P':Pnumber, 'S':Snumber}
    f.close()
    result_list = sorted(result.keys(), key=lambda x: result[x]['adsdelta'])
    return(result, result_list)



def check13C(mz, annotation, threshold = 0.001):
    """Check whether a fragment m/z is the 13C isotope peak of a more intense fragment.

    Scans every other peak in ``annotation`` for one whose m/z is approximately
    one 13C-12C mass unit (``13.00335483534 - 12.0``) below ``mz`` and whose
    intensity is greater than the peak at ``mz`` — i.e. ``mz`` looks like the
    minor 13C isotopologue of that more intense mono-isotopic peak.

    Args:
        mz (float): m/z of the fragment (product ion) peak to check.
        annotation (dict): Mapping of fragment m/z to a dict containing at least
            ``"intensity"``.
        threshold (float, optional): Allowed m/z tolerance around the expected
            13C mass shift. Defaults to ``0.001``.

    Returns:
        bool: ``True`` if ``mz`` is consistent with being the 13C isotope peak of
            a more intense fragment in ``annotation``, ``False`` otherwise.
    """


    targetintensity = float(annotation[mz]["intensity"])

    diff = 13.00335483534 - 12.0

    for fragment_mz in annotation.keys():

        fragmentintensity = float(annotation[fragment_mz]["intensity"])

        if (abs((float(mz) - float(fragment_mz)) - diff) < threshold) and targetintensity < fragmentintensity:
            return True

    return False




def addProton(formula):
    """Add one hydrogen atom to a molecular formula (produce a ``[M+H]`` adduct).

    Args:
        formula (str): Molecular formula of the neutral metabolite.

    Returns:
        str: Molecular formula of the proton-adduct ion.

    See Also:
        :func:`deductProton`: subtraction of one proton.
    """

    mw = ""
    mfs = elementnumber(formula)
    for element, number in mfs.items():
        if element == "H":
            number += 1
        if number == 0:
            continue
        mw = mw+element+str(number)
    return(mw)

def deductProton(formula):
    """Remove one hydrogen atom from a molecular formula (produce a ``[M-H]`` adduct).

    Args:
        formula (str): Molecular formula of the neutral metabolite.

    Returns:
        str: Molecular formula of the proton-deducted ion.

    See Also:
        :func:`addProton`: addition of one proton.
    """

    mf = ""
    mfs = elementnumber(formula)
    for element, number in mfs.items():
        if element == "H":
            number -= 1
        if number <= 0:
            continue
        mf = mf+element+str(number)
    return(mf)

def addAmmonium(formula):
    """Add one ammonium group (NH4) to a molecular formula (produce a ``[M+NH4]`` adduct).

    Args:
        formula (str): Molecular formula of the neutral metabolite.

    Returns:
        str: Molecular formula of the ammonium-adduct ion.

    See Also:
        :func:`addProton`: addition of one proton.
    """
    mw = ""
    mfs = elementnumber(formula)
    for element, number in mfs.items():
        if element == "H":
            number += 4
        if element == "N":
            number += 1
        if number == 0:
            continue
        mw = mw+element+str(number)
    return(mw)

def cleanmf(formula):
    """Normalize a molecular formula string by dropping redundant "1" atom counts.

    For example, ``"H4N1"`` becomes ``"H4N"``. Also removes any elements with a
    zero atom count.

    Args:
        formula (str): Molecular formula to clean up.

    Returns:
        str: Normalized molecular formula.
    """
    mf = ""
    mfs = elementnumber(formula)
    for element, number in mfs.items():
        if number == 0:
            continue
        if number == 1:
            number = ""
        mf = mf+element+str(number)
    return(mf)


def diffformula(mf1, mf2):
    """Compute the elementwise difference of two molecular formulas (``mf1 - mf2``).

    Args:
        mf1 (str): Minuend molecular formula.
        mf2 (str): Subtrahend molecular formula.

    Returns:
        str: Molecular formula representing ``mf1 - mf2``, elementwise. Elements
            whose count difference is exactly zero are omitted (so a formula
            equal to ``mf1 == mf2`` yields ``""``); negative counts are included
            with their sign, e.g. ``"H-2"``.

    See Also:
        :func:`addformula`: addition of two formulas.
    """
    mw = ""
    mfs1 = elementnumber(mf1)
    mfs2 = elementnumber(mf2)
    for element in elements:
        number1 = mfs1[element]
        number2 = mfs2[element]
        if number1-number2 == 0:
            continue

        mw = mw+element+str(number1-number2)
    return(mw)


def addformula(mf1, mf2):
    """Compute the elementwise sum of two molecular formulas (``mf1 + mf2``).

    Args:
        mf1 (str): First molecular formula.
        mf2 (str): Second molecular formula.

    Returns:
        str: Molecular formula representing ``mf1 + mf2``, elementwise. Elements
            whose combined count is exactly zero are omitted.

    See Also:
        :func:`diffformula`: elementwise difference of two formulas.
    """
    mw = ""
    mfs1 = elementnumber(mf1)
    mfs2 = elementnumber(mf2)
    for element in elements:
        number1 = mfs1[element]
        number2 = mfs2[element]
        if number1+number2 == 0:
            continue

        mw = mw+element+str(number1+number2)
    return(mw)

def elementnumber(mf1):
    """Parse a molecular formula into per-element atom counts.

    For every known element in :data:`elements` (Cl, Br, Si, C, H, N, O, S, P, F),
    finds all occurrences of that element symbol in ``mf1`` (each optionally
    followed by a count, defaulting to 1 when omitted) and sums them.

    Args:
        mf1 (str): Molecular formula to parse (e.g. ``"C6H12O6"``).

    Returns:
        dict: Mapping of element symbol (str, from :data:`elements`) to its total
            atom count (int) in ``mf1``; elements absent from the formula map to
            ``0``.
    """
    mf = {}

    for element in elements:
        formula = re.sub(element+r'([^0-9])', element+r'1\1' ,mf1)
        hit = re.findall(element+r'([0-9]*)', mf1)
        newhit = []
        for x in hit:
            if x == '':
                newhit.append(1)
            else:
                newhit.append(x)
        number1 = sum([int(x) for x in newhit])
        mf[element] = number1
    return(mf)


def calcmw(formula, polarity = "none"):
    """Compute the exact (monoisotopic) mass of a molecular formula.

    Sums each element's exact mass (:data:`elements`) times its atom count (via
    :func:`elementnumber`). When ``polarity`` is given, applies the electron-mass
    correction appropriate for a charged ion rather than the neutral molecule.

    Args:
        formula (str): Molecular formula to compute the mass of.
        polarity (str, optional): ``"POSITIVE"`` subtracts one electron mass
            (``0.0005486``); ``"NEGATIVE"`` adds one electron mass; any other
            value (default ``"none"``) returns the neutral-molecule mass
            unmodified.

    Returns:
        float: Exact mass (or exact ion mass, if ``polarity`` is given) of
            ``formula``, in Da.
    """
    mw = 0
    mfs = elementnumber(formula)
    for element, number in mfs.items():
        mw = mw + number * elements[element]

    if polarity == "POSITIVE":
        mw = mw - 0.0005486
    elif polarity == "NEGATIVE":
        mw = mw + 0.0005486
    else:
        mw = mw
    return(mw)

def mzcorrection(mz, correction):
    """Apply a ppm m/z-calibration correction to a measured m/z value.

    ``correction`` may be either a single constant ppm offset, or a string
    encoding a list of ``(mz, correction_ppm)`` calibration points (sorted by
    ascending m/z), in which case the correction level for ``mz`` is obtained by
    piecewise-linear interpolation between the two bracketing calibration points
    when ``mz`` falls strictly between two of them.

    For every other case — ``mz`` below the first point, at or beyond the
    last point, or exactly equal to any calibration point's own m/z — the
    interpolation loop never reaches its ``break`` statement, so the
    ``for...else`` clause always runs afterwards and *overwrites*
    ``correction_level`` with the value linearly extrapolated from just the
    first two calibration points. In particular, the last calibration
    point's own correction value is effectively unreachable: it is only
    ever assigned, never used, because no ``break`` follows that assignment.

    Args:
        mz (float): Measured m/z value to correct.
        correction (float or str): A constant ppm offset, or a string
            representation of a list of ``(mz, correction_ppm)`` tuples
            (evaluated with ``eval``).

    Returns:
        float: Corrected m/z value, computed as
            ``mz + mz * correction_level / 1e6``.
    """
    correction_level = 0.0


    if type(eval(str(correction))) == list:
        correction = eval(correction)
        for i, (cmz1, correct1) in enumerate(correction):
            if cmz1 < mz:
                if (cmz1, correct1) == correction[-1]: # if this is a last
                    correction_level = correct1
                else:
                    (cmz2, correct2) = correction[i+1]
                    if mz < cmz2:
                        correction_level = correct1 + (correct2 - correct1) * (mz - cmz1)/(cmz2 - cmz1)
                        #print(cmz1, correct1, cmz2, correct2, mz, correction_level)
                        break

        else:
            if len(correction) >= 2:
                (cmz1, correct1) = correction[0]
                (cmz2, correct2) = correction[1]
                correction_level = correct1 + (correct1 - correct2) /(cmz1 - cmz2) * (cmz1 - mz)
            else:
                correction_level = 0.0



    else:
        correction_level = float(correction)

    mz_corrected = float(mz) + float(mz) * correction_level /1000000
    #print(mz, mz_corrected, correction_level)
    return(mz_corrected)

def ionform(formula, ionform):
    """Apply a simple adduct/deduct ion-form string to a molecular formula.

    Reads the ion-form's second character as the polarity symbol (``"+"`` for
    addition, anything else for subtraction) and adds or subtracts the remainder
    of the ion-form string (from the 3rd character onward, treated as a formula
    fragment) via :func:`addformula`/:func:`diffformula`.

    Args:
        formula (str): Molecular formula of the neutral metabolite.
        ionform (str): Ion-form string whose 2nd character is the polarity symbol
            (``"+"`` or ``"-"``) and whose remainder (from index 2) is the adduct
            fragment's molecular formula, e.g. ``"[+Na]"``.

    Returns:
        str: Molecular formula of the resulting adduct/deduct ion.

    See Also:
        :func:`generateprecmf`: full-featured precursor-formula generator that
        understands named adducts (``ACN``, ``FA``, ...) and multiply-charged /
        multi-adduct ion-form strings.
    """
    print(ionform[2:])
    pol = ionform[1]
    if pol == "+":
        mf2 = addformula(formula, ionform[2:])
    else:
        mf2 = diffformula(formula, ionform[2:])
    return(mf2)

def generateprecmf(mf, ionform):
    """Compute the precursor-ion molecular formula from a neutral formula and ion-form string.

    Parses a MassBank-style ion-form annotation (e.g. ``"[M+H]+"``,
    ``"[M+FA-H]-"``, ``"[M+Na]+"``) into its individual addition/subtraction
    terms, then applies each term to ``mf`` via :func:`addformula`/
    :func:`diffformula` to build the resulting ion's molecular formula. Named
    adducts (``ACN``, ``FA``, ``Hac``, ``IsoProp``, ``DMSO``) are expanded to
    their underlying chemical formula before being applied; terms consisting
    solely of ``Na`` or ``K`` (charge-carrying alkali-metal adducts) are skipped
    rather than added, so the returned formula reflects the covalently modified
    neutral species rather than the full ionic composition.

    Args:
        mf (str): Molecular formula of the neutral (``M``) species.
        ionform (str): MassBank-style ion-form string, e.g. ``"[M+H]+"``,
            ``"[M-H2O+H]+"``, ``"[M+FA-H]-"``.

    Returns:
        str: Molecular formula of the precursor ion after applying all
            addition/subtraction terms in ``ionform`` (excluding the ignored
            alkali-metal charge carriers).
    """

    adducts ={
    "ACN":"CH3CN",
    "FA": "HCO2H",
    "Hac": "C2H4O2",
    "IsoProp":"C3H8O",
    "DMSO":"C2H6OS",
    }
    adducts_ignore = ["Na", "K"]

    #precmf = genprecmf(mf, ionform)
    ionform = ionform.replace(" ", "")
    ionform = ionform.replace("[", "")
    ionform = ionform.replace("]+", "")
    ionform = ionform.replace("]++", "")
    ionform = ionform.replace("]2+", "")
    ionform = ionform.replace("]-", "")
    ionform = ionform.replace("]--", "")
    ionform = ionform.replace("]2-", "")
    ionform = ionform.replace("]", "")
    ionform = ionform.replace("*", "")

    ionform = ionform.replace("*", "")


    ionform = ionform.replace("+" , "*+")
    ionform = ionform.replace("-" , "*-")

    compositions = ionform.split("*")
    precmf = mf



    for composition in compositions:



        if composition == "M":
            continue
        if "-" in composition:
            composition = composition.replace("-" , "")
            if composition in adducts_ignore:
                continue
            if composition in adducts:
                composition = adducts[composition]

            precmf = diffformula(precmf, composition)

        elif "+" in composition:

            composition = composition.replace("+" , "")
            if composition in adducts_ignore:
                continue
            if composition in adducts:
                composition = adducts[composition]

            precmf = addformula(precmf, composition)
        else:
            if composition in adducts_ignore:
                continue
            if composition in adducts:
                composition = adducts[composition]

            precmf = addformula(precmf, composition)


    #print(ionform, compositions, precmf, mf)


    return(precmf)

