load("./lnp.sage")
load("./katsumata.sage")
load("./utils.sage")

class LNP_Kat(LNP):
    def compute_and_print(self, 
                          mlwe_start_lnp=1,
                          mlwe_samples_lnp=2^10, 
                          mlwe_skip_check_lnp=False,
                          mlwe_start_kat=1, 
                          mlwe_samples_kat=2^10, 
                          mlwe_skip_check_kat=False):
        """
        Compute and check all parameters, and print them.

        Parameters
        ----------
        For LNP
        -------
        mlwe_start_lnp : Initial dimension of MLWE problem.
        mlwe_samples_lnp : Maximum number of samples for the MLWE problems.
        mlwe_skip_check_lnp : Don't use the lattice estimator, just set the MLWE dimension to <mlwe_start_lnp>.

        For Katsumata
        -------------
        mlw_start_kat : Initial dimension of MLWE problem.
        mlwe_samples_kat : Maximum number of samples for the MLWE problems.
        mlwe_skip_check_kat : Don't use the lattice estimator, just set the MLWE dimension to <mlwe_start_kat>.
        """
        print_dashes("LNP")
        super().compute_and_print(mlwe_start=mlwe_start_lnp, 
                                  mlwe_samples=mlwe_samples_lnp, 
                                  mlwe_skip_check=mlwe_skip_check_lnp)
        print_dashes()
        print_dashes("Katsumata")
        total_len = self.m1 + self.m2
        bound = sqrt((self.stdev1*sqrt(2*self.m1*self.d))^2 + (self.stdev2*sqrt(2*self.m2*self.d))^2)
        self.kat = Katsumata(d=self.d, logq=self.logq, m=total_len, B=bound)
        self.kat.compute_and_print(mlwe_start=mlwe_start_kat, 
                                   mlwe_samples=mlwe_samples_kat, 
                                   mlwe_skip_check=mlwe_skip_check_kat)
        print_dashes()
        print_dashes("TOTAL PROOF SIZE")
        self.total_size_with_transform = self.total_size + self.kat.total_size
        print(f"Total proof size: {(self.total_size_with_transform/2^13).n(digits=4)} KB")
        print_dashes()

if __name__ == "__main__":
    lnp_kat = LNP_Kat(d=128, logq=64)
    lnp_kat.compute_and_print(mlwe_start_lnp=20,
                              mlwe_start_kat=21)
    
    lnp_kat = LNP_Kat(d=128, logq=86)
    lnp_kat.compute_and_print(mlwe_start_lnp=27,
                              mlwe_start_kat=28)
    
    lnp_kat = LNP_Kat(d=64, logq=64, omega=8, eta=140, len_s=16, k=16)
    lnp_kat.compute_and_print(mlwe_start_lnp=40,
                              mlwe_start_kat=42)
    
    lnp_kat = LNP_Kat(d=64, logq=90, omega=8, eta=140, len_s=16, k=16)
    lnp_kat.compute_and_print(mlwe_start_lnp=56,
                              mlwe_start_kat=58)