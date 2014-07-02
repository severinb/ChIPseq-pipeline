def execute(cf):
    """
    This function/component refines a given weight matrix, then predicts sites for old and refined WM and true sites and on background peaks.
    With this it makes a positive predictive value - sensitivity plot.
    It also makes a logo of the refined WM.
    """    
    ###
    ### Ports and parameters
    odd_path = cf.get_input("odd_file") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    even_path = cf.get_input("even_file")
    old_wm = cf.get_input("WM")
    WM2 = cf.get_input("WM2")
    bgPeaks = cf.get_input("shuffledPeaks")
    interm = cf.get_output("intermediate")
    ROCplot = cf.get_output("sens_ppv")
    logo = cf.get_output("Logo")
    refWM = cf.get_output("refWM")
    qual_file_nonref = cf.get_output("nonref_WM_Quality")
    qual_file_ref = cf.get_output("ref_WM_Quality")
    
    genome = cf.get_parameter("genome", "string")
    minpostwm = cf.get_parameter("minpostwm", "float")
    WMdiff = cf.get_parameter("wmdiff", "float")
    mylogo_path = cf.get_parameter("mylogo_path", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    ###
    ###
    
    os.mkdir(interm)
    
    wmlen2 = len(open(WM2).readlines())-4    
    os.system('touch ' + logo)     #create output files, so that anduril does not complain
    os.system('touch ' + refWM)
    refWM = WM2  #then put WM2 to refWM, so that further processing works.
    wrong = 0    #unimportant variable, can be set to true or false, makes no difference.
    
    
