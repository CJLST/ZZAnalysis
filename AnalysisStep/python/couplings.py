class Couplings(dict):
   # Spin-0 HQQ couplings
   qqcoupl_spin0 = [
      "kappa",
      "kappa_tilde"
   ]
   # Spin-0 HGG couplings (point-like interaction)
   ggcoupl_spin0 = [
      "ghg2",
      "ghg3",
      "ghg4"
   ]
   # Spin-0 HVV couplings
   cqsq_spin0_HVV = [
      "cz_q1sq",
      "cz_q2sq",
      "cz_q12sq"
   ]
   Lambda_qsq_spin0_HVV = [
      "Lambda_z11",
      "Lambda_z21",
      "Lambda_z31",
      "Lambda_z41",
      "Lambda_z12",
      "Lambda_z22",
      "Lambda_z32",
      "Lambda_z42",
      "Lambda_z10",
      "Lambda_z20",
      "Lambda_z30",
      "Lambda_z40"
   ]
   vvcoupl_spin0_HVV = [
      "ghz1",
      "ghz2",
      "ghz3",
      "ghz4",
      "ghzgs2",
      "ghzgs3",
      "ghzgs4",
      "ghgsgs2",
      "ghgsgs3",
      "ghgsgs4",
      "ghz1_prime",
      "ghz1_prime2",
      "ghz1_prime3",
      "ghz1_prime4",
      "ghz1_prime5",
      "ghz2_prime",
      "ghz2_prime2",
      "ghz2_prime3",
      "ghz2_prime4",
      "ghz2_prime5",
      "ghz3_prime",
      "ghz3_prime2",
      "ghz3_prime3",
      "ghz3_prime4",
      "ghz3_prime5",
      "ghz4_prime",
      "ghz4_prime2",
      "ghz4_prime3",
      "ghz4_prime4",
      "ghz4_prime5",
      "ghzgs1_prime2",
      "ghz1_prime6",
      "ghz1_prime7",
      "ghz2_prime6",
      "ghz2_prime7",
      "ghz3_prime6",
      "ghz3_prime7",
      "ghz4_prime6",
      "ghz4_prime7"
   ]
   cqsq_spin0_HWW = [coupling.replace("z", "w") for coupling in cqsq_spin0_HVV]
   Lambda_qsq_spin0_HWW = [coupling.replace("z", "w") for coupling in Lambda_qsq_spin0_HVV]
   vvcoupl_spin0_HWW = [coupling.replace("z", "w") for coupling in vvcoupl_spin0_HVV]
   # Spin-1 HQQ couplings
   qqcoupl_spin1 = [
      "zprime_qq_left",
      "zprime_qq_right"
   ]
   # Spin-1 HVV couplings
   vvcoupl_spin1 = [
      "zprime_zz_1",
      "zprime_zz_2"
   ]
   # Spin-2 HQQ couplings
   qqcoupl_spin2 = [
      "graviton_qq_left",
      "graviton_qq_right"
   ]
   # Spin-2 HGG couplings
   ggcoupl_spin2 = [
      "a%i"%i for i in range(1,6)
   ]
   # Spin-2 HVV couplings
   vvcoupl_spin2 = [
      "b%i"%i for i in range(1,11)
   ]

   def getspin(self, doreweighting):
       if not doreweighting: return 0
       spin = {S: False for S in (0, 1, 2)}
       for coupling in self.qqcoupl_spin0 + self.ggcoupl_spin0 + self.cqsq_spin0_HVV + self.Lambda_qsq_spin0_HVV + self.vvcoupl_spin0_HVV + self.cqsq_spin0_HWW + self.Lambda_qsq_spin0_HWW + self.vvcoupl_spin0_HWW:
           if self[coupling]:
               spin[0] = True
       for coupling in self.qqcoupl_spin1 + self.vvcoupl_spin1:
           if self[coupling]:
               spin[1] = True
       for coupling in self.qqcoupl_spin2 + self.ggcoupl_spin2 + self.vvcoupl_spin2:
           if self[coupling]:
               spin[2] = True

       if sum(spin.values()) > 1:
           raise ValueError("Couplings for more than one spin set (" + ", ".join(k for k, v in spin.iteritems() if v) + ").  Please set couplings for only one spin value: ghz* for spin 0, zprime_zz_* for spin 1, or a* and b* for spin 2.")

       if sum(spin.values()) == 0:
           raise ValueError("Reweighting is turned on, but no couplings are set!")

       for k, v in spin.iteritems():
           if v:
               return k

   def __setitem__(self, key, value):
      try:
         if isinstance(value, str):   #allow writing 3+2i, where python wants 3+2j
            value = value.replace("i", "j")
         if (("Lambda_z" in key) or ("Lambda_w" in key)):
            value = float(value)
         elif (("cz_" in key) or ("cw_" in key)):
            value = int(value)
         else:
            value = complex(value)
      except ValueError:
         raise ValueError("The coupling {} needs to be an integer, a float or a complex number!".format(key))
      super(Couplings, self).__setitem__(key, value)

   def allnames(self):
      return self.qqcoupl_spin0 + self.ggcoupl_spin0 + self.cqsq_spin0_HVV + self.Lambda_qsq_spin0_HVV + self.vvcoupl_spin0_HVV + self.cqsq_spin0_HWW + self.Lambda_qsq_spin0_HWW + self.vvcoupl_spin0_HWW + self.qqcoupl_spin1 + self.vvcoupl_spin1 + self.qqcoupl_spin2 + self.ggcoupl_spin2 + self.vvcoupl_spin2

   def getcouplings(self, key, imag=False):
      if key == "Hqqcoupl":
         if not imag:
            return [self[name].real for name in self.qqcoupl_spin0]
         else:
            return [self[name].imag for name in self.qqcoupl_spin0]
      elif key == "Hggcoupl":
         if not imag:
            return [self[name].real for name in self.ggcoupl_spin0]
         else:
            return [self[name].imag for name in self.ggcoupl_spin0]
      elif key == "HzzCLambda_qsq":
         return [self[name].real for name in self.cqsq_spin0_HVV]
      elif key == "HzzLambda_qsq":
         return [self[name].real for name in self.Lambda_qsq_spin0_HVV]
      elif key == "Hzzcoupl":
         if not imag:
            return [self[name].real for name in self.vvcoupl_spin0_HVV]
         else:
            return [self[name].imag for name in self.vvcoupl_spin0_HVV]
      elif key == "HwwCLambda_qsq":
         return [self[name].real for name in self.cqsq_spin0_HWW]
      elif key == "HwwLambda_qsq":
         return [self[name].real for name in self.Lambda_qsq_spin0_HWW]
      elif key == "Hwwcoupl":
         if not imag:
            return [self[name].real for name in self.vvcoupl_spin0_HWW]
         else:
            return [self[name].imag for name in self.vvcoupl_spin0_HWW]
      elif key == "Zqqcoupl":
         if not imag:
            return [self[name].real for name in self.qqcoupl_spin1]
         else:
            return [self[name].imag for name in self.qqcoupl_spin1]
      elif key == "Zvvcoupl":
         if not imag:
            return [self[name].real for name in self.vvcoupl_spin1]
         else:
            return [self[name].imag for name in self.vvcoupl_spin1]
      elif key == "Gqqcoupl":
         if not imag:
            return [self[name].real for name in self.qqcoupl_spin2]
         else:
            return [self[name].imag for name in self.qqcoupl_spin2]
      elif key == "Gggcoupl":
         if not imag:
            return [self[name].real for name in self.ggcoupl_spin2]
         else:
            return [self[name].imag for name in self.ggcoupl_spin2]
      elif key == "Gvvcoupl":
         if not imag:
            return [self[name].real for name in self.vvcoupl_spin2]
         else:
            return [self[name].imag for name in self.vvcoupl_spin2]
      assert False

