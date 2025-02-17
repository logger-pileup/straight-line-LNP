load("./utils.sage")

class SL_LNP:
    """
    Represents a straight-line LNP proof configuration.

    We compute the non-interactive proof size of proving knowledge of a MLWE sample.
    Concretely, we want to prove knowledge of a vector s such that the norm of s, 
    ||(s,Ds-u)||, is at most B.
    """
    def __init__(self,
                 d = 128, 
                 logq = 86,
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
                 gamma3 = 2):
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
        v : maximum coefficent size of commitment randomness.
        k : height of D.
        B : exact norm bound being proved.
        eta : the heuristic bound on \sqrt[2k](|| \sigma_{-1}(c^k)*c^k ||_1) for k = 32
        gamma1 : rejection sampling for z1 (masked small part of commitment)
        gamma2 : rejection sampling for z2 (masked commitment randomness)
        gamma3 : rejection sampling for Rs^(e) (ARP)
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
        self.gamma3 = gamma3

        # Find an exact value for q.
        # We need q to be prime and q = 2l+1 (mod 4l).
        self.q = 2^(self.logq) + (2*self.l + 1)        
        self.q =  self.q - 4*self.l                      # find the next candidate for q
        while not is_prime(self.q):
            self.q -= 4*self.l

        # number of repetitions for boosting soundness, we assume lambda is even
        self.lmbda = 2*ceil(secpam/(2*logq))

        # length of pBs2 + (s1||m) part of the commitment
        self.m1 = len_s + 1                              # short part of commitment contains s and x := bin. decomp. of B^2 - ||(s,As-u)||^2                                                     # 
        self.ell = 256/d + 1 + self.lmbda/2 + 1          # large part of commitment contains y3 (mask in ARP), a bit (in ARP), lambda/2 garbadge polynomials (in eval) and g1 (in single quad)
        self.mdash = self.m1 + 2*self.ell                # the total length of the pBs2+(s1||m) part of the commitment

        # rejection sampling
        self.alphae = sqrt(2048 + d)                     # bound alpha^(e) on the vector e^(e) = (s,As - u, bin. decomp. of B^2 - ||(s,As-u)||^2)
        self.stdev1 = gamma1 * eta * sqrt(alpha^2 + d)   # s.d. for masking short part of commitment 
        self.stdev3 = gamma3 * sqrt(337) * self.alphae   # s.d. for y3 which masks opening in ARP

        # store any errors
        self.errors = []

        # soundness error
        self.soundness_err = 2/(2*omega+1)^(d/2) + self.lmbda/(self.q^(d/2)) + (258+self.d)/self.q^(self.lmbda) +2*d*(self.len_s+self.k+1)/2^256
        if self.soundness_err > 2^-secpam:
            self.errors.append(f"Parameters do not meet requirements for soundness. Soundness error: 2^{log(self.soundness_err, 2).n()}")

        # flags
        self.mlwe_complete = False                       # has the MLWE dimension been found
        self.requirements_checked = False                # requirements for extraction have been checked
        self.size_set = False                            # have set the size

    def compute_mlwe(self, start=1, samples=2^10, skip_check=False):
        """
        We need to ensure that the MLWE problem is hard to ensure hiding and
        extractability of the commitment.

        For extraction, we have mdash problems in parallel (which gives 
        reduction loss of 1/mdash).

        We find large enough dimensions so that commitment is both hiding and
        extractable. We set the security to be the worse of the two cases.

        The lattice estimator can be slow, so the skip_check parameter can be
        used to set the MLWE dimension to <start>.

        Parameters
        ----------
        start : Initial dimension of MLWE problem.
        samples : Maximum number of samples for the MLWE problems.
        skip_check : Don't use the lattice estimator, just set the MLWE dimension to <start>.
        """
        if self.mlwe_complete:
            print(f"Already set the MLWE dimensions: hid={self.mlwe_dim_hid}, ext={self.mlwe_dim_ext}")
            return
        
        mlwe_dim = start-1                                  # dimension being increased
        self.mlwe_samples = samples                         # number of samples
        self.mlwe_dim_hid = 0                               # dimension for hiding problem, set when found
        self.mlwe_dim_ext = 0                               # dimension for extraction problem, set when found
        self.mlwe = 0                                       # MLWE security

        print("Finding the MLWE dimension...")
        while self.mlwe < self.secpam:
            mlwe_dim += 1

            if skip_check:
                print(f"Skipping lattice estimator for LWE - setting the dimension to: {start}")
                self.mlwe_dim_hid = self.mlwe_dim_ext = mlwe_dim
                break

            mlwe = lwe_estimate(
                n=mlwe_dim, d=self.d, logq=self.logq, nu=self.v, m=self.mlwe_samples
            )
            print(f"MLWE estimation results: dim={mlwe_dim}: {mlwe} bits")

            # if this dimension is good enough for hiding, then set it
            if mlwe >= self.secpam and not self.mlwe_dim_hid:
                self.mlwe_dim_hid = mlwe_dim

            # adjust for reduction loss for repeating the problem m times in parallel
            mlwe_parallel = mlwe - ceil(log(self.mdash, 2))

            # if this dimension is good enough for extraction, then set it
            if mlwe_parallel >= self.secpam and not self.mlwe_dim_ext:
                self.mlwe_dim_ext = mlwe_dim
            
            # set the secuirty as the minimum of the two
            self.mlwe = min(mlwe, mlwe_parallel)

        # check that the the number of samples has not been exceeded
        if self.mlwe_dim_hid + self.mlwe_dim_ext + self.mdash > self.mlwe_samples:
            self.errors.append(f"Cannot guarantee security of MLWE for hiding since the commitment key width \
                                exceeds the specified number of samples")

        if self.mlwe_dim_ext + self.mdash > self.mlwe_samples:
            self.errors.append(f"Cannot guarantee security of MLWE for hiding since the commitment \
                                length exceeds the specified number of samples")

        # set the height and width of A
        self.m = self.mlwe_dim_ext
        self.m2 = self.mlwe_dim_hid + self.mlwe_dim_ext + self.mdash

        # set s.d. for rejection sampling on commitment randomness
        self.stdev2 = self.gamma2 * self.eta * self.v * sqrt(self.m2 * self.d)

        # repitition rate
        t = sqrt(2*(self.secpam+1)/log(e,2))
        self.rep_rate = exp(t/self.gamma1 + 1/(2*self.gamma1^2)) * exp(t/self.gamma2 + 1/(2*self.gamma2^2)) * exp(1/(2*self.gamma3^2))

        self.mlwe_complete = True
    
    def check_requirements(self):
        """
        Check the requirements on extraction and on ARP.

        We require that p and q are large enough so that no wrap occurs
        during extraction with z' and r' (Lemma 3.1). This ensures the 
        binding (unique extraction) of the witness.

        We also need to ARP proof to ensure no wrap occurs mod q in the
        later relations.
        """
        if self.requirements_checked:
            print("Requirements for extraction already checked.")
            return
        
        if not self.mlwe_complete:
            print("You must set the MLWE dimension ")
            return

        print("Setting p and checking requirements...")
        
        self.norm_z1_prime = 8*self.d*self.omega*self.stdev1*sqrt(2*self.m1)       # infinity norm of z1 prime (Lemma 3.1)
        self.norm_z2_prime = 8*self.d*self.omega*self.stdev2*sqrt(2*self.m2)       # infinity norm of z2 prime (Lemma 3.1)

        # Set value of p
        self.p = ceil(2*self.norm_z1_prime+1)

        # Estimate how many bits q needs to be
        estimate_q = max(2*(self.p*self.v*self.norm_z2_prime*self.m2*self.d + self.norm_z1_prime), 
                         2*self.v*self.norm_z2_prime*self.m2*self.d*(1+sqrt(self.q)))
        print(f"Estimate that q should be at least {log(estimate_q, 2).n(digits=4)} bits")

        # Check requirements for extraction of z prime
        if self.p*self.v*self.norm_z2_prime*self.m2*self.d + self.norm_z1_prime >= self.q/2:
            self.errors.append("Cannot extract small part of z prime - wraps mod q")

        if self.norm_z1_prime >= self.p/2:
            self.errors.append("Cannot extract small part of z prime - wraps mod p")

        if self.v*self.norm_z2_prime*self.m2*self.d*(1+sqrt(self.q)) >= self.q/2:
            self.errors.append("Cannot extract large part of z prime - wraps mod q")

        if self.v*self.norm_z2_prime*self.m2*self.d >= sqrt(self.q)/2:
            self.errors.append("Cannot extract large part of z prime - wraps mod root q")

        # Check ARP ensures no wrap mod q
        arp_norm = 2 * sqrt((self.len_s+self.k+1) * self.d) * 1.64 * self.stdev3

        if self.q <= arp_norm^2 + arp_norm*sqrt(self.d):
            self.errors.append("Cannot prove all x_i have binary coefficients - ARP does not guarantee no wrap mod q")

        if self.q <= 3 * self.B^2 + arp_norm^2:
            print("Cannot prove ||s,Ds-u|| <= B - ARP does not guarantee no wrap mod q")

        self.requirements_checked = True

    def calculate_size(self):
        """
        Calculate the proof size.

        The proof size consists of
            - the commitment size
            - size of all other full-sized polynomials that need to be sent (lambda/2 of them)
            - size of all small polynomials (masked openings)
            - size of the challenge
        """
        if self.size_set:
            print("Have already calculated the size")
            return
        
        if not self.requirements_checked:
            print("Cannot calculate the size without calculating the MLWE dimension and checking requirements first")
            return
    
        self.com_size = (self.m + self.mdash) * self.d * self.logq                      # total commitment size
        self.full_pols = self.lmbda/2 * self.d * self.logq                              # masked opening of garbage polynomials for eval
        self.short_pols = self.m2 * self.d * ceil(log(self.stdev2,2) + 2.57)            # masked opening of commitment randomness
        self.short_pols += self.m1 * self.d * (ceil(log(self.stdev1,2) + 2.57))         # masked opening of m1 := s || bin. decomp. of B^2 - ||(s,As-u)||^2
        self.short_pols += 256 * (ceil(log(self.stdev3, 2) + 2.57))                     # masked opening of y3 (ARP)
        self.chal_size = ceil(log(2*self.omega+1,2)) * self.d                           # size of the challenge
        self.total_size = self.com_size + self.full_pols + self.short_pols + self.chal_size
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
        print(f"Rejection sampling: gamma1={self.gamma1}, gamma2={self.gamma2}, gammae={self.gamma3}")
        print_dashes()

    def print_computed_parameters(self):
        """
        Print the paramaters that have been computed.
        """
        print_dashes("COMPUTED PARAMETERS")
        print(f"Exact modulus: {self.q}")

        if self.requirements_checked:
            print(f"Commitment factor p: {self.p} ({log(self.p, 2).n(digits=4)} bits)")
        
        print(f"Lambda for boosting soundness: {self.lmbda}")
        print(f"Length of short part of committed elements: {self.m1}")
        print(f"Length of long part of committed elements: {self.ell}")
        print(f"Height of A: {self.m if self.mlwe_complete else 'not set'}")
        print(f"Height of B: {self.mdash}")
        print(f"Length of commitment randomness: {self.m2 if self.mlwe_complete else 'not set'}")
        print(f"Rejection sampling: stdev1={self.stdev1.n(digits=8)}, stdev2={self.stdev2.n(digits=8) if self.mlwe_complete else 'not set'}, stdev3={self.stdev3.n(digits=8)}")
        print_dashes()

    def print_security_analysis(self):
        """
        Print security analysis.
        """
        print_dashes("SECURITY ANALYSIS")
        print(f"Soundness Error: 2^{(ceil(log(self.soundness_err, 2).n()))}")

        if self.mlwe_complete:
            print(f"MLWE dimension for extraction: {self.mlwe_dim_ext}")
            print(f"MLWE dimension for hiding: {self.mlwe_dim_hid}")
            print(f"MLWE security: {self.mlwe} bits")
            print(f"MLWE samples: {self.mlwe_samples}")
            print(f"Repitition rate: {log(self.rep_rate, 2).n(digits=4)} bits")

        if self.requirements_checked:
            print(f"Infinity norm z1 prime: {self.norm_z1_prime.n(digits=8)} ({log(self.norm_z1_prime, 2).n(digits=4)}) bits")
            print(f"Infinity norm z2 prime: {self.norm_z2_prime.n(digits=8)} ({log(self.norm_z2_prime, 2).n(digits=4)}) bits")

        print_dashes()

    def print_size(self):
        """
        Print the size.
        """
        print_dashes("SIZE")

        if self.size_set:
            print(f"Commitment size: {(self.com_size/2^13).n(digits=4)} KB")
            print(f"Short polynomials size: {(self.short_pols/2^13).n(digits=4)} KB")
            print(f"Full-sized polynomials size: {(self.full_pols/2^13).n(digits=4)} KB")
            print(f"Total size: {(self.total_size/2^13).n(digits=4)} KB")
            
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
        
        if not self.requirements_checked:
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
        mlw_start : Initial dimension of MLWE problem.
        mlwe_samples : Maximum number of samples for the MLWE problems.
        mlwe_skip_check : Don't use the lattice estimator, just set the MLWE dimension to <start>.
        """
        self.compute_mlwe(mlwe_start, mlwe_samples, mlwe_skip_check)
        self.check_requirements()
        self.calculate_size()
        self.print()

if __name__ == "__main__":
    sl_lnp = SL_LNP(d=128, logq=86, omega=2, eta=59)
    sl_lnp.compute_and_print(mlwe_start=27)

    sl_lnp = SL_LNP(d=64, logq=90, omega=8, eta=140, len_s=16, k=16)
    sl_lnp.compute_and_print(mlwe_start=56)
