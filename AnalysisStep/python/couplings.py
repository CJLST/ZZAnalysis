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
    def __setitem__(self, key, value):
        try:
            if isinstance(value, str):   #allow writing 3+2i, where python wants 3+2j
                value = value.replace("i", "j")
            value = complex(value)
        except ValueError:
            raise ValueError("coupling {} needs to be a complex number!".format(key))
        super(Couplings, self).__setitem__(key, value)

    def getcouplings(self, WW=False, imag=False):
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
