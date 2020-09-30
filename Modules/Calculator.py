from __future__ import print_function
# -*- coding: utf-8 -*-
import cellconstructor as CC
import cellconstructor.Methods



import difflib
__ASE_CALC__ = False
try:
    import ase
    from ase.calculators.espresso import Espresso
    __ASE_CALC__ = True
except:
    pass

"""
This class setup the ASE calculator for the run from
a standard namelist.
"""

__CALCULATOR_HEAD__ = "calculator"
__KPTS_HEAD__ = "k_points"
__KOFF_HEAD__ = "k_offset"
__DISABLE_CHECK__ = "disable_check"
__CALCULATOR_TYPE__ = "program"
__BINARY__ = "binary"
__PSEUDO_LIST__ = "pseudo_"

__CALC_QE__  = "qe"
__CALCULATOR_SYNONIMOUS__ = {"qe" : __CALC_QE__, 
                             "quantumespresso" : __CALC_QE__,
                             "quantum-espresso" : __CALC_QE__}

__QE_ALLOWED_KEYS__ = ["ecutrho", "ecutwfc", "smearing", "degauss", 
                       "occupations", "conv_thr", "tstress", "tprnfor",
                       "verbosity", "disk_io", "input_dft", "use_all_frac"]

__ALLOWED_KEYS__ = [__KOFF_HEAD__, __KPTS_HEAD__, __DISABLE_CHECK__, 
                    __CALCULATOR_TYPE__, __BINARY__]

__REQUESTED_KEYS__ = [__CALCULATOR_TYPE__, __KPTS_HEAD__]

def prepare_calculator_from_namelist(namelist):
    """
    PREPARE ASE CALCULATOR
    ======================
    
    This function prepares the ASE calculator for the execution.
    It requires a namelist for the setup.
    
    Parameters
    ----------
        namelist : dict or string
            The parsed namelist.
            If a string is passed, it will be parsed
    
    Returns
    -------
        ase_calc :
            The ASE calculator.
    """

    # Parse the namelist if needed
    if isinstance(namelist, str):
        namelist = CC.Methods.read_namelist(namelist)
    
    # Check if the namelist has the correct keys
    if not __CALCULATOR_HEAD__ in namelist.keys():
        raise ValueError("Error, to setup a calculator the section %s must be declared." % __CALCULATOR_HEAD__)
    
    c_info = namelist[__CALCULATOR_HEAD__]
    ks = c_info.keys()
    
    # Check if the disable check is present
    check = True
    if __DISABLE_CHECK__ in ks:
        check = bool(c_info[__DISABLE_CHECK__])
    
    # Define defaults
    KPTS = (1,1,1)
    KOFFS = (0,0,0)
    
    
    # Get the keywords
    if __KPTS_HEAD__ in ks:
        KPTS = [int(x) for x in c_info[__KPTS_HEAD__]]
        if len(KPTS) != 3:
            raise ValueError("Error, %s must be a 3 dimensional list" % __KPTS_HEAD__)
    
    if __KOFF_HEAD__ in ks:
        KOFFS = [int(x) for x in c_info[__KOFF_HEAD__]]
        if len(KOFFS) != 3:
            raise ValueError("Error, %s must be a 3 dimensional list" % __KOFF_HEAD__)
    
    # Setup the pseudopotentials
    pseudopotentials = {}
    pseudo_keys = [x for x in ks if __PSEUDO_LIST__ in x]
    for pkey in pseudo_keys:
        pvalue = c_info[pkey]
        patom = pkey.replace(__PSEUDO_LIST__, "", 1)
        patom = patom[0].upper() + patom[1:]
        pseudopotentials[patom] = pvalue
        
    is_binary = False
    if __BINARY__ in ks:
        is_binary = True
        binary = c_info[__BINARY__]
    
    # Setup the calculator
    ase_calc = None
    tot_keys = __ALLOWED_KEYS__
    if __CALCULATOR_TYPE__ in ks:
        if not c_info[__CALCULATOR_TYPE__] in __CALCULATOR_SYNONIMOUS__.keys():
            print ("List of supported calculators:", __CALCULATOR_SYNONIMOUS__.keys)
            raise ValueError("Error, the specified calculator '%s' is not in supported." % c_info[__CALCULATOR_TYPE__])
            
        calc = __CALCULATOR_SYNONIMOUS__[c_info[__CALCULATOR_TYPE__]]
        
        
        if calc == __CALC_QE__:            
            tot_keys = __ALLOWED_KEYS__ + __QE_ALLOWED_KEYS__
            # Check all the qe allowed keys
            input_data = {}
            qe_keys = [x for x in __QE_ALLOWED_KEYS__ if x in ks]
            for x in qe_keys:
                input_data[x] = c_info[x]
                
            ase_calc = Espresso(pseudopotentials = pseudopotentials,
                                input_data = input_data, kpts = KPTS,
                                koffset = KOFFS)
            
            if is_binary:
                ase_calc.command = binary
    
    # Check the allowed keys
    for k in ks: 
        if __PSEUDO_LIST__ in k:
            continue
        if (not k in tot_keys) and check:
            print ("Error with the key:", k)
            print ("Did you mean something like:", difflib.get_close_matches(k, tot_keys))
            raise IOError("Error in calculator namespace: key '" + k +"' not recognized.")
    
    # Check for mandatory keys
    for req_key in __REQUESTED_KEYS__:
        if not req_key in ks:
            raise IOError("Error, the calculator configuration namelist requires the keyword: '" + req_key + "'")
         
    
    return ase_calc
    
