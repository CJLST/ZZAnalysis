float update_xsec(TString inputname, float xsec) {
        float new_xsec;

        // Updated BR values.
        // xsec values same as our csv files, as they are correct.
        // We use xsec*BR to weight the events, hence here we recompute it.

        if (inputname.Contains("bbH120")) { new_xsec = 0.5534*0.0001659; }
        else if (inputname.Contains("bbH124")) { new_xsec = 0.4999*0.0002502; }
        else if (inputname.Contains("bbH125")) { new_xsec = 0.4880*0.0002745; }
        else if (inputname.Contains("bbH126")) { new_xsec = 0.4760*0.0003001; }
        else if (inputname.Contains("bbH130")) { new_xsec = 0.4304*0.0004124; }

        else if (inputname.Contains("ttH120")) { new_xsec = 0.5697*0.00042015426186342315; }
        else if (inputname.Contains("ttH124")) { new_xsec = 0.5193*0.000649581819587205; }
        else if (inputname.Contains("ttH125") || inputname.Contains("ttH125_tuneup") || inputname.Contains("ttH125_tunedown")) { new_xsec = 0.5071*0.0007176792246182684; }
        else if (inputname.Contains("ttH126")) { new_xsec = 0.4964*0.000788372579311756; }
        else if (inputname.Contains("ttH130")) { new_xsec = 0.4539*0.0010947747168867862; }

        else if (inputname.Contains("ggZH125")) { new_xsec = 0.1227*0.0006991275831543255; }

        else if (inputname.Contains("ZH120")) { new_xsec = 0.993920889*0.0004106877051958817;}
        else if (inputname.Contains("ZH124")) { new_xsec = 0.90514981*0.0006322344128856341;}
        else if (inputname.Contains("ZH125") || inputname.Contains("ZH125_tuneup") || inputname.Contains("ZH125_tunedown")) { new_xsec = 0.7612*0.0006991275831543255;}
        else if (inputname.Contains("ZH126")) { new_xsec = 0.865060831*0.0007650556125370585;}
        else if (inputname.Contains("ZH130")) { new_xsec = 0.789702354*0.0010674878938974918;}

        else new_xsec = xsec;

        return new_xsec;
}
