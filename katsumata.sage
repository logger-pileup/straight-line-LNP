load("./utils.sage")

class Katsumata:
    """
    Represents a Katsumata transform protocol.

    The Katsumata transform bootstraps many lattice-based ZKPs
    into straight-line extractable ZKPs.

    In order to extract multiple transcripts without rewinding, we
    commit to the witness with an extractable, linear-homomorphic 
    (LinHC) commitment scheme, and additionally send a masked opening
    of the LinHC commitment randomness. This can then be used to
    extract more masked openings of the witness without rewinding.
    """
    def __init__(self, 
                 d = 128, 
                 logq = 64,
                 secpam = 128,
                 m=8,
                 B=0,
                 omega = 2, 
                 v = 1, 
                 gamma = 14):
        """
        Construct on the configurable parameters.

        Parameters
        ----------
        d : dimension on Rq.
        logq : number of bits of modulus, q.
        secpam : security parameter.
        m : length of the vector that needs to be committed to.
        B: L2 norm bound on vector that needs to be committed to.
        omega : maximum coefficient size of challenges.
        v : maximum coefficent size of commitment randomness and of trapdoor matrices.
        gamma: rejection sampling parameter for commitment randomness.
        """
        # Set configured parameters
        self.d = d
        self.logq = logq
        self.secpam = secpam
        self.m = m
        self.B = B
        self.omega = omega
        self.v = v
        self.gamma = gamma

        # store any errors
        self.errors = []

        # flags
        self.mlwe_complete = False          # has the MLWE dimension been found
        self.requirements_checked = False   # have the requirements been checked
        self.size_set = False               # have set the size

    def compute_mlwe(self, start=1, samples=2^10, skip_check=False):
        """
        We need to ensure that the MLWE problem is hard to ensure hiding and
        extractability of the commitment.

        For hiding, we have dimension n. For extraction, we have dimension m and
        m cases in parallel (which gives reduction loss of 1/m).

        We find the dimension n which is large enough, taking into account the 
        reduction loss, and then check n <= m.

        The lattice estimator can be slow, so the skip_check parameter can be
        used to set the MLWE dimension to <start>.

        Parameters
        ----------
        start : Initial dimension of MLWE problem.
        samples : Maximum number of samples for the MLWE problems.
        skip_check : Don't use the lattice estimator, just set the MLWE dimension to <start>.
        """
        if self.mlwe_complete:
            print(f"Already set the MLWE dimension: {self.n}")
            return
        
        self.n = start-1
        self.mlwe_samples = samples
        self.mlwe = 0

        print("Finding the MLWE dimension...")
        while self.mlwe < self.secpam:
            self.n += 1

            if skip_check:
                print(f"Skipping lattice estimator for LWE - setting the dimension to: {start}")
                break

            self.mlwe = lwe_estimate(
                n=self.n, d=self.d, logq=self.logq, nu=self.v, m=self.mlwe_samples
            )
            self.mlwe -= ceil(log(self.m, 2))                                               # reduction loss for repeating the problem m times in parallel
            print(f"MLWE estimation results: dim={self.n}: {self.mlwe} bits")

        if self.n > self.m:
            self.errors.append(f"Cannot gurantee security of extraction - MLWE dimension n ({self.n}) exceeds length m ({self.m})")

        # check that the the number of samples has not been exceeded
        if self.m > self.mlwe_samples or self.n > self.mlwe_samples:
            self.errors.append(f"Cannot guarantee security of MLWE for hiding since the specified number of samples may have been exceeded")

        # rejection sampling
        self.stdev = self.gamma * self.v * self.omega * sqrt((self.n + 2*self.m) * self.d)   # s.d. for masking commitment randomness

        self.mlwe_complete = True
    
    def check_requirements(self):
        """
        Check the requirements extraction.

        We require q to be set such that the norm of the computed element from extraction
        does not wrap mod q.

        We require p to be greater than twice the norm of the element being extracted so
        that it does not wrap mod p.
        """
        if self.requirements_checked:
            print("Requirements for extraction already checked.")
            return
        
        if not self.mlwe_complete:
            print("You must set the MLWE dimension ")
            return

        print("Setting p and checking requirements...")

        # Set value of p
        self.p = ceil(4*self.B+1)
        if not is_odd(self.p):
            self.p += 1

        # Estimate how many bits q needs to be
        norm = 2*sqrt(2)*self.p*(self.n*self.d*self.v + sqrt(self.n*self.m) * self.d * self.v + sqrt(self.n*self.d)) * self.stdev + 2*self.B
        print(f"Estimate that q should be at least {log(2*norm, 2).n(digits=4)} bits")

        # Check requirements for extraction
        if self.B < sqrt(2*self.n*self.d) * self.stdev:
            self.errors.append("Norm bound on committed vector is too small")

        if norm >= 2^(self.logq-1):
            self.errors.append("Cannot extract - may wrap mod q")

        if self.B > (self.p-1)/4:
            self.errors.append("Cannot extract - may wrap mod p")

        self.requirements_checked = True

    def calculate_size(self):
        """
        Calculate the size of the transform in bits.

        The size will be the sum of the commitment, masking opening of the
        commitment randomness and challenge.
        """
        if self.size_set:
            print("Have already calculated the size")
            return
        
        if not self.requirements_checked:
            print("Cannot calculate the size without calculating the MLWE dimension and checking requirements first")
            return
        
        self.com_size = 2 * self.m * self.d * self.logq 
        self.op_size = (self.n + 2*self.m) * self.d * log(10 * self.stdev, 2)
        self.chal_size = ceil(log(2*self.omega+1,2)) * self.d
        self.total_size = self.com_size + self.op_size + self.chal_size
        self.size_set = True
        
    def print_input_parameters(self):
        """
        Print the parameters that were provided.
        """
        print_dashes("INPUT PARAMETERS")
        print(f"Ring dimension: {self.d}")
        print(f"Number of bits of modulus: {self.logq}")
        print(f"Security parameter: {self.secpam}")
        print(f"Length of committed vector: {self.m}")
        print(f"Norm bound on committed vector: {self.B.n(digits=8)} ({log(self.B, 2).n(digits=4)} bits)")
        print(f"Maximum coefficient size of challenge: {self.omega}")
        print(f"Maximum coefficient size of commitment randomness: {self.v}")
        print(f"Rejection sampling: gamma={self.gamma}")
        print_dashes()

    def print_computed_parameters(self):
        """
        Print the paramaters that have been computed.
        """
        print_dashes("COMPUTED PARAMETERS")
        if self.mlwe_complete:
            print(f"n: {self.n}")
            print(f"Rejection sampling: stdev={self.stdev.n(digits=8)}")

        if self.requirements_checked:
            print(f"p: {self.p.n(digits=8)} ({log(self.p, 2).n(digits=4)} bits)")
    
        print_dashes()

    def print_security_analysis(self):
        """
        Print security analysis.
        """
        print_dashes("SECURITY ANALYSIS")

        if self.mlwe_complete:
            print(f"MLWE dimension: {self.n}")
            print(f"MLWE security: {self.mlwe} bits")
            print(f"MLWE samples: {self.mlwe_samples}")

        print_dashes()

    def print_size(self):
        """
        Print the size.
        """
        print_dashes("SIZE")

        if self.size_set:
            print(f"Commitment size: {(self.com_size/2^13).n(digits=4)} KB")
            print(f"Open size: {(self.op_size/2^13).n(digits=4)} KB")
            print(f"Challenge: {(self.chal_size/2^13).n(digits=4)} KB")
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
