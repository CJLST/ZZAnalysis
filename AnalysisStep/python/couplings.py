class Couplings(dict):
    couplingsorder_Z = [
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
                        "ghz4_prime7",
                       ]
    couplingsorder_W = [coupling.replace("z", "w") for coupling in couplingsorder_Z]
    couplingsorder_spin1 = [
                            "zprime_zz_1",
                            "zprime_zz_2",
                           ]
    couplingsorder_spin2_gg = [
                               "a%i"%i for i in range(1,6)
                              ]
    couplingsorder_spin2 = [
                            "b%i"%i for i in range(1,11)
                           ]
    def __setitem__(self, key, value):
        try:
            if isinstance(value, str):   #allow writing 3+2i, where python wants 3+2j
                value = value.replace("i", "j")
            value = complex(value)
        except ValueError:
            raise ValueError("coupling {} needs to be a complex number!".format(key))
        super(Couplings, self).__setitem__(key, value)

    def allnames(self):
        return self.couplingsorder_Z + self.couplingsorder_W + self.couplingsorder_spin1 + self.couplingsorder_spin2_gg + self.couplingsorder_spin2

    def getspin(self, doreweighting):
        if not doreweighting: return 0
        spin = {S: False for S in (0, 1, 2)}
        for coupling in self.couplingsorder_Z + self.couplingsorder_W:
            if self[coupling]:
                spin[0] = True
        for coupling in self.couplingsorder_spin1:
            if self[coupling]:
                spin[1] = True
        for coupling in self.couplingsorder_spin2_gg + self.couplingsorder_spin2:
            if self[coupling]:
                spin[2] = True

        if sum(spin.values()) > 1:
            raise ValueError("Couplings for more than one spin set (" + ", ".join(k for k, v in spin.iteritems() if v) + ").  Please set couplings for only one spin value: ghz* for spin 0, zprime_zz_* for spin 1, or a* and b* for spin 2.")

        if sum(spin.values()) == 0:
            raise ValueError("Reweighting is turned on, but no couplings are set!")

        for k, v in spin.iteritems():
            if v:
                return k

    def getcouplings(self, spin, gg=False, WW=False, imag=False):
        if spin == 0:
            assert not gg
            if not WW:
                if not imag:
                    return [self[name].real for name in self.couplingsorder_Z]
                else:
                    return [self[name].imag for name in self.couplingsorder_Z]
            else:
                if not imag:
                    return [self[name].real for name in self.couplingsorder_W]
                else:
                    return [self[name].imag for name in self.couplingsorder_W]

        if spin == 1:
            assert not gg and not WW
            if not imag:
                return [self[name].real for name in self.couplingsorder_spin1]
            else:
                return [self[name].imag for name in self.couplingsorder_spin1]

        if spin == 2:
            assert not WW
            if not gg:
                if not imag:
                    return [self[name].real for name in self.couplingsorder_spin2]
                else:
                    return [self[name].imag for name in self.couplingsorder_spin2]
            else:
                if not imag:
                    return [self[name].real for name in self.couplingsorder_spin2_gg]
                else:
                    return [self[name].imag for name in self.couplingsorder_spin2_gg]

        assert False
