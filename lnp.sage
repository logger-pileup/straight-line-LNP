load("./utils.sage")

class LNP:
    """
    Represents a LNP proof configuration.

    We compute the non-interactive proof size of proving knowledge of a MLWE sample.
    Concretely, we want to prove knowledge of a vector s such that the norm of s, 
    ||(s,Ds-u)||, is at most B.
    """
    def __init__(self, 
                 d = 128, 
                 logq = 32,
                 secpam = 128,
                 l = 2,
                 omega=2, 
                 len_s=8, 
                 alpha = sqrt(1024), 
                 v = 1,
                 k = 8, 
                 B = sqrt(2048),
                 eta = 59,
                 gamma1 = 19,
                 gamma2 = 2,
                 gammae = 2):
        """
        Construct on the configurable parameters.

        Parameters
        ----------
        d : dimension on Rq.
        logq : number of bits of modulus, q.
        secpam : security parameter.
        l : number of irreducible factors of X^d + 1 modulo q.
        omega : maximum coefficient size of challenges.
        len_s : length of witness s.
        alpha : L2 norm bound on s.
        v : maximum coefficent size of ABDLOP commitment randomness.
        k : height of D.
        B : exact norm bound being proved.
        eta : the heuristic bound on \sqrt[2k](|| \sigma_{-1}(c^k)*c^k ||_1) for k = 32
        gamma1 : rejection sampling for z1 (masked small part of commitment)
        gamma2 : rejection sampling for z2 (masked commitment randomness)
        gammae : rejection sampling for Rs^(e) (ARP)
        """
        # Set configured parameters
        self.d = d
        self.logq = logq
        self.secpam = secpam
        self.l = l
        self.omega = omega
        self.len_s = len_s
        self.alpha = alpha
        self.v = v
        self.k = k
        self.B = B
        self.eta = eta
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gammae = gammae

        # number of repetitions for boosting soundness, we assume lambda is even
        self.lmbda = 2*ceil(secpam/(2*logq))

        # length of commitment
        self.m1 = len_s + 1                              # short part of commitment contains s and x := bin. decomp. of B^2 - ||(s,As-u)||^2
        self.ell = self.lmbda/2 + 256/self.d + 2         # long part of commitment contains lambda/2 (eval) + 256/d + 1 (ARP) + 1 (quad-single) garbage polynomials
        
        # rejection sampling
        self.alphae = sqrt(2048 + d)                     # bound alpha^(e) on the vector e^(e) = (s,As - u, bin. decomp. of B^2 - ||(s,As-u)||^2)
        self.stdev1 = gamma1 * eta * sqrt(alpha^2 + d)   # s.d. for masking short part of commitment 
        self.stdeve = gammae * sqrt(337) * self.alphae   # s.d. for y3 which masks opening in ARP

        # store any errors
        self.errors = []

        # flags
        self.mlwe_complete = False                       # has the MLWE dimension been found
        self.msis_complete = False                       # has the MSIS dimension been found
        self.requirements_checked = False                # have the requirements been checked
        self.size_set = False                            # have set the size
        self.error = False                               # no satisfactory parameters can be found

    def compute_mlwe(self, start=11, samples=2^10, skip_check=False):
        """
        We need to ensure that MLWE is hard to ensure the commitment is hiding.

        Parameters
        ----------
        start : Initial dimension of MLWE problem.
        samples : Maximum number of samples for the MLWE problems.
        skip_check : Don't use the lattice estimator, just set the MLWE dimension to <start>.
        """
        if self.error:
            return
        
        if self.mlwe_complete:
            print(f"Already set the MLWE dimension: {self.mlwe_dim}")
            return

        self.mlwe_dim = start-1                         # dimension of problem
        self.mlwe_samples = samples                     # mlwe samples
        self.mlwe = 0                                   # mlwe security

        print("Finding MLWE dimension...")
        while self.mlwe < self.secpam:                          
            self.mlwe_dim += 1                          # increasing the mlwe dimension until MLWE provides ~ 128-bit security

            if skip_check:
                print(f"Skipping lattice estimator for LWE - setting the dimension to: {self.mlwe_dim}")
                break
            
            self.mlwe = lwe_estimate(
                n=self.mlwe_dim, d=self.d, logq=self.logq, nu=self.v, m=self.mlwe_samples
            )
            print(f"MLWE estimation results: mlwe-dim={self.mlwe_dim}: {self.mlwe} bits")

        self.mlwe_complete = True

    def compute_msis(self):
        """
        We need to esnure that MSIS is hard to ensure the commitment is bidning.

        Additionally, we find the parameters for the commitment compression techniques
        from Dilithium-G [Duc+17].
        """
        if self.error:
            return

        if self.msis_complete:
            print(f"Already set the MSIS dimension: {self.msis_dim}")
            return

        if not self.mlwe_complete:
            print("Please compute the MLWE dimension before the MSIS dimension.")
            return

        print("Finding MSIS dimension...")
        self.msis_dim = 0                                                                   # dimension of the Module-SIS problem
        self.D = 0                                                                          # dropping low-order bits of t_A
        self.gamma = 0                                                                      # dropping low-order bits of w
        self.msis = 0                                                                       # hardness of MSIS
        self.msis_bound = 0                                                                 # bound on the extracted MSIS solution
        Bound1 =  2 * self.stdev1 * sqrt(2 * (self.m1 + 1) * self.d)                        # bound on bar{z}_1

        # Find the SIS dimension to meet the required security                                                                       
        while self.msis < self.secpam:                                                            
            self.msis_dim += 1                                                              # inrease dimension until required security is met

            self.m2 = self.mlwe_dim + self.msis_dim + self.ell                              # length of commitment randomness           
            self.stdev2 = self.gamma2 * self.eta * self.v * sqrt(self.m2 * self.d)          # set stdev2 with the current candidate for msis_dim
    
            Bound2 =  (2 * self.stdev2 * sqrt(2 * self.m2 * self.d) + 2^self.D * self.eta   # bound on bar{z}_2 = (bar{z}_{2,1},bar{z}_{2,2})
                       * sqrt(self.msis_dim*self.d) + self.gamma * sqrt(self.msis_dim*self.d))   
            self.msis_bound = 4 * self.eta * sqrt(Bound1^2 + Bound2^2)                      # bound on the extracted MSIS solution

            # The bound can't exceed the modulus. If it does, we need a bigger modulus.
            if self.msis_bound >= (2^self.logq - 1)/2:
                self.errors.append("The bound on the extracted MSIS solution exceeds the modulus. No parameters can be found, use a bigger modulus.")
                self.error = True
                return

            self.msis = sis_estimate(self.msis_dim, self.d, self.logq, self.msis_bound)     # estimate SIS hardness

        # check that the the number of samples has not been exceeded
        if self.m2 - self.mlwe_dim > self.mlwe_samples:
            self.errors.append("Number of MLWE samples has been exceeded by the size of the commitment")
        
        # Given the SIS dimension, find the largest possible gamma which makes the MSIS solution still small
        print("Computing the parameter gamma...")
        self.msis = 0
        self.gamma = 2^self.logq                                                            # initialisation
        while self.msis < self.secpam:                                                      # searching for right gamma
            self.gamma /= 2                                                                 # decrease the value of gamma    
            Bound2 =  (2 * self.stdev2 * sqrt(2 * self.m2 * self.d) + 2^self.D * self.eta   # bound on bar{z}_2 = (bar{z}_{2,1},bar{z}_{2,2})
                       * sqrt(self.msis_dim*self.d) + self.gamma * sqrt(self.msis_dim*self.d))   
            self.msis_bound = 4 * self.eta * sqrt(Bound1^2 + Bound2^2)                      # bound on the extracted MSIS solution

            # The bound can't exceed the modulus
            if self.msis_bound >= (2^self.logq - 1)/2:
                continue

            self.msis = sis_estimate(self.msis_dim, self.d, self.logq, self.msis_bound)     # estimate SIS hardness

        # Now we can find the exact value of q and gamma
        # We need q to be prime and q = 2l+1 (mod 4l).
        # We need gamma to be a divisor of q-1
        print("Computing exact q and gamma...")
        true_gamma_found = False                                                            # Boolean for finding correct gamma
        self.q = 2^(self.logq) + (2*self.l + 1)                                             # we need q to be congruent to 2l+1 modulo 4l
        while not true_gamma_found:
            self.q -= 4*self.l                                                              # find the next candidate for q
            while not is_prime(self.q):                                                     # we need q to be prime 
                self.q -= 4*self.l

            div_q = divisors(self.q-1)                                                      # consider divisors of q-1
            for div in div_q:                 
                if self.gamma*4/5 < div and div <= self.gamma and is_even(div):             # find a divisor which is close to gamma
                    self.gamma = div                                                        # we found a good candidate for gamma
                    true_gamma_found = True

        # soundness error
        self.soundness_err = 2 * 1/(2*self.omega+1)^(self.d/2) +  self.q^(-self.d/self.l) + self.q^(-self.lmbda) + 2^(-128) 
        if self.soundness_err > 2^-self.secpam:
            self.errors.append(f"Parameters do not meet requirements for soundness. Soundness error: 2^{log(self.soundness_err, 2).n()}")

        # Given msis dimension and gamma, find the largest possible D which makes the MSIS solution small
        print("Computing the parameter D...")
        self.msis = 0
        self.D = self.logq
        while self.msis < self.secpam: 
            self.D -= 1                                                                     # decrease the value of D    
            Bound1 =  2 * self.stdev1 * sqrt(2 * (self.m1 + 1) * self.d)                    # bound on bar{z}_1
            Bound2 =  (2 * self.stdev2 * sqrt(2 * self.m2 * self.d) + 2^self.D * self.eta   # bound on bar{z}_2 = (bar{z}_{2,1},bar{z}_{2,2})
                       * sqrt(self.msis_dim*self.d) + self.gamma * sqrt(self.msis_dim*self.d))   
            self.msis_bound = 4 * self.eta * sqrt(Bound1^2 + Bound2^2)                      # bound on the extracted MSIS solution
            
            if self.msis_bound >= (2^self.logq - 1)/2 or 2^(self.D-1) * self.omega * self.d >= self.gamma:
                continue

            self.msis = sis_estimate(self.msis_dim, self.d, self.logq, self.msis_bound)     # estimate SIS hardness

        # repitition rate
        t = sqrt(2*(self.secpam+1)/log(e,2))
        self.rep_rate = 2*exp(t/self.gamma1 + 1/(2*self.gamma1^2)) * exp(1/(2*self.gamma2^2)) * exp(1/(2*self.gammae^2))
        self.msis_complete = True
    
    def check_requirements(self):
        """
        Check knowledge soundness conditions from Theorem 5.3 (LNP).
        """
        if self.error:
            return
        
        if self.requirements_checked:
            print("Requirements for extraction already checked.")
            return
        
        if not self.msis_complete:
            print("You must set the MLWE and MSIS dimensions first")
            return
        
        print("Checking knowledge soundness conditions...")

        Be = 2 * sqrt(256/26) * 1.64 * self.stdeve
        if self.q <  41 * (self.len_s + self.k + 1) * self.d * Be:
            self.errors.append("Cannot use Lemma 2.9")

        if self.q <= Be^2 + Be*sqrt(self.d):
            self.errors.append("Cannot prove all x_i have binary coefficients")

        if self.q <= 3 * self.B^2 + Be^2:
            self.errors.append("Cannot prove || E_i*s - v_i || <= beta_i")

        self.requirements_checked = True

    def calculate_size(self):
        """
        Calculate the proof size.

        The proof size consists of
            - the commitment size
            - size of all other full-sized polynomials that need to be sent (lambda/2 of them)
            - size of all small polynomials (masked openings)
            - size of the challenge
            - size of the hint (for commitment compression)
        """
        if self.error:
            return
        
        if self.size_set:
            print("Have already calculated the size")
            return
        
        if not self.requirements_checked:
            print("Cannot calculate the size without calculating the MLWE dimension and checking requirements first")
            return

        self.com_size = self.msis_dim * self.d * (self.logq - self.D) + self.ell * self.d * self.logq   # total commitment size (using compression for short part)
        self.full_pols = self.lmbda/2 * self.d * self.logq                                              # masked opening of garbage polynomials for eval
        self.short_pols = 256 * (ceil(log(self.stdeve, 2) + 2.57))                                      # masked opening of y3 (ARP)
        self.short_pols += (self.m2-self.msis_dim) * self.d * ceil(log(self.stdev2,2) + 2.57)           # masked opening of m2
        self.short_pols += self.m1 * self.d * (ceil(log(self.stdev1,2) + 2.57))                         # masked opening of m1 := s || bin. decomp. of B^2 - ||(s,As-u)||^2
        self.chal_size = ceil(log(2*self.omega+1,2)) * self.d                                           # size of the challenge
        self.hint_size = 2.25 * self.msis_dim * self.d
        self.total_size = self.com_size + self.full_pols + self.short_pols + self.chal_size + self.hint_size
        self.size_set = True
        
    def print_input_parameters(self):
        """
        Print the parameters that were provided.
        """
        print_dashes("INPUT PARAMETERS")
        print(f"Ring dimension: {self.d}")
        print(f"Number of bits of modulus: {self.logq}")
        print(f"Security parameter: {self.secpam}")
        print(f"Number of irreducible factors of Rq: {self.l}")
        print(f"Maximum coefficient size of challenge: {self.omega}")
        print(f"Length of witness s: {self.len_s}")
        print(f"Norm bound on s: {self.alpha}")
        print(f"Maximum coefficient size of commitment randomness: {self.v}")
        print(f"Height of D: {self.k}")
        print(f"Exact norm bound to prove: {self.B.n(digits=4)}")
        print(f"Challenge set: eta={self.eta}")
        print(f"Rejection sampling: gamma1={self.gamma1}, gamma2={self.gamma2}, gammae={self.gammae}")
        print_dashes()

    def print_computed_parameters(self):
        """
        Print the paramaters that have been computed.
        """
        if self.error:
            return
        
        print_dashes("COMPUTED PARAMETERS")       
        print(f"Exact modulus: {self.q}")        
        print(f"Lambda for boosting soundness: {self.lmbda}")
        print(f"Length of short part of committed elements: {self.m1}")
        print(f"Length of long part of committed elements: {self.ell}")

        if self.msis_complete:
            print(f"Length of commitment randomness: {self.m2}")
            print(f"Parameter gamma for dropping low-order bits of w: {self.gamma}")
            print(f"Parameter D for dropping low-order bits of t_A: {self.D}")
            print(f"Rejection sampling: stdev1={self.stdev1.n(digits=8)}, stdev2={self.stdev2.n(digits=8)}, stdeve={self.stdeve.n(digits=8)}")
        print_dashes()

    def print_security_analysis(self):
        """
        Print security analysis.
        """
        if self.error:
            return
        
        print_dashes("SECURITY ANALYSIS")
        print(f"Soundness Error: 2^{(log(self.soundness_err, 2).n(digits=4))}")

        if self.mlwe_complete:
            print(f"MLWE dimension: {self.mlwe_dim}")
            print(f"MLWE security: {self.mlwe.n(digits=4)} bits")
            print(f"MLWE sample: {self.mlwe_samples}")

        if self.msis_complete:
            print(f"SIS dimension: {self.msis_dim}")
            print(f"SIS security: {self.msis.n(digits=4)} bits")
            print(f"Bound on the extracted MSIS solution: {self.msis_bound.n(digits=8)}")
            print(f"Repitition rate: {self.rep_rate.n(digits=4)} bits")

        print_dashes()

    def print_size(self):
        """
        Print the proof size.
        """
        if self.error:
            return
        
        print_dashes("PROOF SIZE")

        if self.size_set:
            print(f"Commitment: {(self.com_size/2^13).n(digits=4)} KB")
            print(f"Other full-sized polynomials: {(self.full_pols/2^13).n(digits=4)} KB")
            print(f"Short polynomials: {(self.short_pols/2^13).n(digits=4)} KB")
            print(f"Challenge: {(self.chal_size/2^13).n(digits=4)} KB")
            print(f"Hint: {(self.hint_size/2^13).n(digits=4)} KB")
            print(f"Total proof size: {(self.total_size/2^13).n(digits=4)} KB")

        print_dashes()

    def print_errors(self):
        """
        Print errors.
        """
        print_dashes("ERRORS")

        for err in self.errors:
            print(err)

        if not self.mlwe_complete:
            print("WARNING: Have not found MLWE dimensions")
        
        if not self.error and not self.requirements_checked:
            print("WARNING: Have not checked requirements")
        print_dashes()

    def print(self):
        """
        Print all.
        """
        self.print_input_parameters()
        self.print_computed_parameters()
        self.print_security_analysis()
        self.print_size()
        self.print_errors()

    def compute_and_print(self, mlwe_start=1, mlwe_samples=2^10, mlwe_skip_check=False):
        """
        Compute and check all parameters, and print them.

        Parameters
        ----------
        mlwe_start : Initial dimension of MLWE problem.
        mlwe_samples : Maximum number of samples for the MLWE problems.
        mlwe_skip_check : Don't use the lattice estimator, just set the MLWE dimension to <mlwe_start>.
        """
        self.compute_mlwe(mlwe_start, mlwe_samples, mlwe_skip_check)
        self.compute_msis()
        self.check_requirements()
        self.calculate_size()
        self.print()

if __name__ == "__main__":
    print("If you intended to run LNP+Katsumata, please remove this block")
    lnp = LNP(d=128, logq=32)
    lnp.compute_and_print(mlwe_start=11)

    lnp = LNP(d=128, logq=64)
    lnp.compute_and_print(mlwe_start=20)

    lnp = LNP(d=128, logq=86)
    lnp.compute_and_print(mlwe_start=27)

    lnp = LNP(d=64, logq=64, omega=8, eta=140, len_s=16, k=16)
    lnp.compute_and_print(mlwe_start=40)

    lnp = LNP(d=64, logq=90, omega=8, eta=140, len_s=16, k=16)
    lnp.compute_and_print(mlwe_start=56)
