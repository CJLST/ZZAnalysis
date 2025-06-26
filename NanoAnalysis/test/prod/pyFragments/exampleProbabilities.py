 setConf("probabilities", {'Name': "Native", 
                          "Process": "SelfDefine_spin0", 
                          "MatrixElement":"JHUGen",
                          "Production": "ZZGG", 
                          "Couplings": {'ghz1':[1,0], 
                                        'ghg2':[1,0]}, 
                          "Prod": False,
                          "Dec": True, 
                          "isgen": True, 
                          "computeprop": False 
                          }, 
                          append=True)

    setConf("probabilities", {'Name': "P_Gen_ggH_ghg2_1_ghz4_1", 
                            "Process": "SelfDefine_spin0", 
                            "MatrixElement":"JHUGen",
                            "Production": "ZZGG", 
                            "Couplings": {'ghz4':[1,0], 
                                            'ghg2':[1,0]}, 
                            "Prod": False,
                            "Dec": True, 
                            "isgen": True,
                            "dividep": "Native",
                            "computeprop": False }, append=True)