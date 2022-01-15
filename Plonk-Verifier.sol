pragma solidity >=0.5.0 <0.7.0;

library PairingsBn254 {
    uint256 constant q_mod = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 constant r_mod = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    
    struct G1Point {
        uint256 X;
        uint256 Y;
    } 
    
    struct Fr {
        uint256 value;
    }
    
    function new_fr(uint256 fr) internal pure returns (Fr memory) {
        require(fr < r_mod);
        return Fr({value: fr});
    }
    
    function copy(Fr memory self) internal pure returns (Fr memory n) {
        n.value = self.value;
    }
    
    function assign(Fr memory self, Fr memory other) internal pure {
        self.value = other.value;
    }
    
    function inverse(Fr memory fr) internal view returns (Fr memory) {
        assert(fr.value != 0);
        return pow(fr, r_mod-2);
    }
    
    function add_assign(Fr memory self, Fr memory other) internal pure {
        self.value = addmod(self.value, other.value, r_mod);
    }
    
    function sub_assign(Fr memory self, Fr memory other) internal pure {
        self.value = addmod(self.value, r_mod - other.value, r_mod);
    }
    
    function mul_assign(Fr memory self, Fr memory other) internal pure {
        self.value = mulmod(self.value, other.value, r_mod);
    }
    
    function pow(Fr memory self, uint256 power) internal view returns (Fr memory) {
        uint256[6] memory input = [32, 32, 32, self.value, power, r_mod];
        uint256[1] memory result;
        bool success;
        assembly {
            success := staticcall(gas(), 0x05, input, 0xc0, result, 0x20)
        }
        require(success);
        return Fr({value: result[0]});
    }
    
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }

    function P1() internal pure returns (G1Point memory) {
        return G1Point(1, 2);
    }
    
    function new_g1(uint256 x, uint256 y) internal pure returns (G1Point memory) {
        return G1Point(x, y);
    }
    
    function new_g2(uint256[2] memory x, uint256[2] memory y) internal pure returns (G2Point memory) {
        return G2Point(x, y);
    }
    
    function copy_g1(G1Point memory self) internal pure returns (G1Point memory result) {
        result.X = self.X;
        result.Y = self.Y;
    }

    function P2() internal pure returns (G2Point memory) {
        // for some reason ethereum expects to have c1*v + c0 form
        
        return G2Point(
            [0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2,
                0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed],
            [0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b,
                0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa]
        );
    }

    function negate(G1Point memory self) internal pure {
        // The prime q in the base field F_q for G1
        if (self.X == 0 && self.Y == 0)
            return;
        self.Y = q_mod - self.Y;
    }
    
    // function is_infinity(G1Point memory p) internal pure returns (bool) {
    //     if (p.X == 0 && p.Y == 0) {
    //         return true;
    //     }
    //     return false;
    // }

    function point_add(G1Point memory p1, G1Point memory p2)
        internal view returns (G1Point memory r)
    {
        point_add_into_dest(p1, p2, r);
        return r;
    }
    
    function point_add_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        point_add_into_dest(p1, p2, p1);
    }

    function point_add_into_dest(G1Point memory p1, G1Point memory p2, G1Point memory dest)
        internal view
    {
        uint256[4] memory input;
        if (p2.X == 0 && p2.Y == 0) {
            // we add zero, nothing happens
            dest.X = p1.X;
            dest.Y = p1.Y;
            return;
        } else if (p1.X == 0 && p1.Y == 0) {
            // we add into zero, and we add non-zero point
            dest.X = p2.X;
            dest.Y = p2.Y;
            return;
        } else {
            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = p2.Y;
        }
        bool success = false;
        assembly {
            success := staticcall(gas(), 6, input, 0x80, dest, 0x40)
        }
        require(success);
    }
    
    function point_sub_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        point_sub_into_dest(p1, p2, p1);
    }

    function point_sub_into_dest(G1Point memory p1, G1Point memory p2, G1Point memory dest)
        internal view
    {
        uint256[4] memory input;
        if (p2.X == 0 && p2.Y == 0) {
            // we subtracted zero, nothing happens
            dest.X = p1.X;
            dest.Y = p1.Y;
            return;
        } else if (p1.X == 0 && p1.Y == 0) {
            // we subtract from zero, and we subtract non-zero point
            dest.X = p2.X;
            dest.Y = q_mod - p2.Y;
            return;
        } else {
            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = q_mod - p2.Y;
        }
        bool success = false;
        assembly {
            success := staticcall(gas(), 6, input, 0x80, dest, 0x40)
        }
        require(success);
    }


    function point_mul(G1Point memory p, Fr memory s)
        internal view returns (G1Point memory r)
    {
        point_mul_into_dest(p, s, r);
        return r;
    }
    
    function point_mul_assign(G1Point memory p, Fr memory s)
        internal view
    {
        point_mul_into_dest(p, s, p);
    }

    function point_mul_into_dest(G1Point memory p, Fr memory s, G1Point memory dest)
        internal view
    {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s.value;
        bool success;
        assembly {
            success := staticcall(gas(), 7, input, 0x60, dest, 0x40)
        }
        require(success);
    }
    
    function pairing(G1Point[] memory p1, G2Point[] memory p2)
        internal view returns (bool)
    {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := staticcall(gas(), 8, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
        }
        require(success);
        return out[0] != 0;
    }

    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2)
        internal view returns (bool)
    {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
}

library TranscriptLibrary {
    // flip                    0xe000000000000000000000000000000000000000000000000000000000000000;
    uint256 constant FR_MASK = 0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff;

    uint32 constant DST_0 = 0;
    uint32 constant DST_1 = 1;
    uint32 constant DST_CHALLENGE = 2;
    
    struct Transcript {
        bytes hash_data;
        uint32 challenge_counter;
    }

    function new_transcript() internal pure returns (Transcript memory t) {
        t.hash_data = abi.encodePacked();
        t.challenge_counter = 0;
    }

    function update_with_u256(Transcript memory self, uint256 value) internal pure {
        self.hash_data = abi.encodePacked(self.hash_data, value);
    }
    
    function update_with_fr(Transcript memory self, PairingsBn254.Fr memory value) internal pure {
        update_with_u256(self, value.value);
    }
    
    function update_with_g1(Transcript memory self, PairingsBn254.G1Point memory p) internal pure {
        update_with_u256(self, p.X);
        update_with_u256(self, p.Y);
    }
    
    function get_challenge(Transcript memory self, uint256 second) internal pure returns(PairingsBn254.Fr memory challenge) {
        bytes32 query = sha256(abi.encodePacked(self.hash_data));
        if(second != 0) {
            query = sha256(abi.encodePacked(self.hash_data, uint256(query) & FR_MASK));
        }
        self.challenge_counter += 1;
        challenge = PairingsBn254.Fr({value: uint256(query) & FR_MASK});
    }
}

contract Plonk4VerifierWithAccessToDNext {
    using PairingsBn254 for PairingsBn254.G1Point;
    using PairingsBn254 for PairingsBn254.G2Point;
    using PairingsBn254 for PairingsBn254.Fr;
    
    using TranscriptLibrary for TranscriptLibrary.Transcript;

    uint256 constant STATE_WIDTH = 3;
    
    struct VerificationKey { // "Verifier preprocessed input"
        uint256 domain_size; // n
        uint256 num_inputs; // len(public input)
        PairingsBn254.Fr omega; // omega
        // STATE_WIDTH for witness + multiplication + constant
        PairingsBn254.G1Point[STATE_WIDTH+2] selector_commitments; // [q_L.R.O.M.C]_1
        PairingsBn254.G1Point[STATE_WIDTH] permutation_commitments; // [s_sigma1.2.3]_1
        PairingsBn254.Fr[STATE_WIDTH-1] permutation_non_residues; // k1 k2
        PairingsBn254.G2Point g2_x; // [x]_2
    }
    
    struct Proof {
        uint256[] input_values; // public inputs w[]
        PairingsBn254.G1Point[STATE_WIDTH] wire_commitments; // [a.b.c]_1
        PairingsBn254.G1Point grand_product_commitment; // [z]_1
        PairingsBn254.G1Point[STATE_WIDTH] quotient_poly_commitments; // [t_lo.mid.hi]_1
        PairingsBn254.Fr[STATE_WIDTH] wire_values_at_z; // _a.b.c
        PairingsBn254.Fr grand_product_at_z_omega; // _z_omega
        PairingsBn254.Fr quotient_polynomial_at_z; // _t; in std-plonk this is computed by V
        PairingsBn254.Fr linearization_polynomial_at_z; // _r
        PairingsBn254.Fr[STATE_WIDTH-1] permutation_polynomials_at_z; // _s_sigma1.2
    
        PairingsBn254.G1Point opening_at_z_proof; // [W_z]_1
        PairingsBn254.G1Point opening_at_z_omega_proof; // [W_z_omega]_1
    }
    
    struct PartialVerifierState {
        PairingsBn254.Fr alpha;
        PairingsBn254.Fr beta;
        PairingsBn254.Fr gamma;
        PairingsBn254.Fr v;
        PairingsBn254.Fr u;
        PairingsBn254.Fr z;
        PairingsBn254.Fr[] cached_lagrange_evals; // L_i(z)
    }
    
    function evaluate_lagrange_poly_out_of_domain( // L_i(z); i = polynum
        uint256 poly_num, 
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        require(poly_num < domain_size);
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);
        PairingsBn254.Fr memory omega_power = omega.pow(poly_num);
        res = at.pow(domain_size);
        res.sub_assign(one);
        assert(res.value != 0); // Vanishing polynomial can not be zero at point `at`
        res.mul_assign(omega_power);
        
        PairingsBn254.Fr memory den = PairingsBn254.copy(at);
        den.sub_assign(omega_power);
        den.mul_assign(PairingsBn254.new_fr(domain_size));
        
        den = den.inverse();
        res.mul_assign(den);
    }
    
    function likely_batch_evaluate_lagrange_poly_out_of_domain( // matter-labs' batch function has problem
		uint256[] memory poly_nums, // 0 ~ n-1
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
	) internal view returns (PairingsBn254.Fr[] memory res) {
		PairingsBn254.Fr[] memory nums = new PairingsBn254.Fr[](poly_nums.length);
		for (uint i = 0; i < poly_nums.length; i++) {
            nums[i] = evaluate_lagrange_poly_out_of_domain(poly_nums[i], domain_size, omega, at);
        }
        return nums;
	}
    
    function evaluate_vanishing( // "zero polynomial evaluation" zeta^n-1
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        res = at.pow(domain_size);
        res.sub_assign(PairingsBn254.new_fr(1));
    }
    
    function verify_at_z( // this function check _t because it is provided by P now; should we modify to fit std-plonk and zkrypt output?
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        PairingsBn254.Fr memory lhs = evaluate_vanishing(vk.domain_size, state.z);
        assert(lhs.value != 0); // we can not check a polynomial relationship if point `z` is in the domain
        lhs.mul_assign(proof.quotient_polynomial_at_z);
    
        PairingsBn254.Fr memory quotient_challenge = PairingsBn254.new_fr(1);
        PairingsBn254.Fr memory rhs = PairingsBn254.copy(proof.linearization_polynomial_at_z);
        
        // public inputs
        PairingsBn254.Fr memory tmp = PairingsBn254.new_fr(0);
        for (uint256 i = 0; i < proof.input_values.length; i++) {
            tmp.assign(state.cached_lagrange_evals[i]);
            tmp.mul_assign(PairingsBn254.new_fr(proof.input_values[i]));
            rhs.add_assign(tmp);
        }
        
        quotient_challenge.mul_assign(state.alpha);
        
        PairingsBn254.Fr memory z_part = PairingsBn254.copy(proof.grand_product_at_z_omega);

        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp.assign(proof.permutation_polynomials_at_z[i]);
            tmp.mul_assign(state.beta);
            tmp.add_assign(state.gamma);
            tmp.add_assign(proof.wire_values_at_z[i]);
            
            z_part.mul_assign(tmp);
        }
        
        tmp.assign(state.gamma);
        // we need a wire value of the last polynomial in enumeration
        tmp.add_assign(proof.wire_values_at_z[STATE_WIDTH - 1]);
        
        z_part.mul_assign(tmp);
        z_part.mul_assign(quotient_challenge);
        
        rhs.sub_assign(z_part);
        
        quotient_challenge.mul_assign(state.alpha);
        
        tmp.assign(state.cached_lagrange_evals[0]); // L_1(z)
        tmp.mul_assign(quotient_challenge);
        
        rhs.sub_assign(tmp);
        
        return lhs.value == rhs.value;
    }
    
    function reconstruct_d( // compute [D]_1
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (PairingsBn254.G1Point memory res) {
    	
        res = PairingsBn254.copy_g1(vk.selector_commitments[STATE_WIDTH + 1]);
                
        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();
        PairingsBn254.Fr memory tmp_fr = PairingsBn254.new_fr(0);
        
        // addition gates
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            tmp_g1 = vk.selector_commitments[i].point_mul(proof.wire_values_at_z[i]);
            res.point_add_assign(tmp_g1);
        }
        
        // multiplication gate
        tmp_fr.assign(proof.wire_values_at_z[0]);
        tmp_fr.mul_assign(proof.wire_values_at_z[1]);
        tmp_g1 = vk.selector_commitments[STATE_WIDTH].point_mul(tmp_fr);
        res.point_add_assign(tmp_g1);
        
        // z * non_res * beta + gamma + a
        PairingsBn254.Fr memory grand_product_part_at_z = PairingsBn254.copy(state.z);
        grand_product_part_at_z.mul_assign(state.beta);
        grand_product_part_at_z.add_assign(proof.wire_values_at_z[0]);
        grand_product_part_at_z.add_assign(state.gamma);
        for (uint256 i = 0; i < vk.permutation_non_residues.length; i++) {
            tmp_fr.assign(state.z);
            tmp_fr.mul_assign(vk.permutation_non_residues[i]);
            tmp_fr.mul_assign(state.beta);
            tmp_fr.add_assign(state.gamma);
            tmp_fr.add_assign(proof.wire_values_at_z[i+1]);
            
            grand_product_part_at_z.mul_assign(tmp_fr);
        }
        
        grand_product_part_at_z.mul_assign(state.alpha);
    
        tmp_fr.assign(state.cached_lagrange_evals[0]);
        tmp_fr.mul_assign(state.alpha);
        tmp_fr.mul_assign(state.alpha);
        
        grand_product_part_at_z.add_assign(tmp_fr);
        
        PairingsBn254.Fr memory last_permutation_part_at_z = PairingsBn254.new_fr(1);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp_fr.assign(state.beta);
            tmp_fr.mul_assign(proof.permutation_polynomials_at_z[i]);
            tmp_fr.add_assign(state.gamma);
            tmp_fr.add_assign(proof.wire_values_at_z[i]);
            
            last_permutation_part_at_z.mul_assign(tmp_fr);
        }

        last_permutation_part_at_z.mul_assign(state.beta);
        last_permutation_part_at_z.mul_assign(proof.grand_product_at_z_omega);
        last_permutation_part_at_z.mul_assign(state.alpha);
        
        // add to the linearization
        tmp_g1 = proof.grand_product_commitment.point_mul(grand_product_part_at_z);
        tmp_g1.point_sub_assign(vk.permutation_commitments[STATE_WIDTH - 1].point_mul(last_permutation_part_at_z));

        res.point_add_assign(tmp_g1);
        res.point_mul_assign(state.v);
        
        res.point_add_assign(proof.grand_product_commitment.point_mul(state.u));
    }
    
    function verify_commitments(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        PairingsBn254.G1Point memory d = reconstruct_d(state, proof, vk);
        
        PairingsBn254.Fr memory z_in_domain_size = state.z.pow(vk.domain_size);
        
        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();

        PairingsBn254.Fr memory aggregation_challenge = PairingsBn254.new_fr(1);
        
        // compute [F]_1 as "commitment_aggregation"
        PairingsBn254.G1Point memory commitment_aggregation = PairingsBn254.copy_g1(proof.quotient_poly_commitments[0]);
        PairingsBn254.Fr memory tmp_fr = PairingsBn254.new_fr(1);
        for (uint i = 1; i < proof.quotient_poly_commitments.length; i++) {
            tmp_fr.mul_assign(z_in_domain_size);
            tmp_g1 = proof.quotient_poly_commitments[i].point_mul(tmp_fr);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        aggregation_challenge.mul_assign(state.v);
        commitment_aggregation.point_add_assign(d);
        
        for (uint i = 0; i < proof.wire_commitments.length; i++) {
            aggregation_challenge.mul_assign(state.v);
            tmp_g1 = proof.wire_commitments[i].point_mul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }
        
        for (uint i = 0; i < vk.permutation_commitments.length - 1; i++) {
            aggregation_challenge.mul_assign(state.v);
            tmp_g1 = vk.permutation_commitments[i].point_mul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }
        
        // collect opening values and compute [E]_1 as "aggregation_challenge"
        aggregation_challenge = PairingsBn254.new_fr(1);
        
        PairingsBn254.Fr memory aggregated_value = PairingsBn254.copy(proof.quotient_polynomial_at_z);
        
        aggregation_challenge.mul_assign(state.v);

        tmp_fr.assign(proof.linearization_polynomial_at_z);
        tmp_fr.mul_assign(aggregation_challenge);
        aggregated_value.add_assign(tmp_fr);
        
        for (uint i = 0; i < proof.wire_values_at_z.length; i++) {
            aggregation_challenge.mul_assign(state.v);
            
            tmp_fr.assign(proof.wire_values_at_z[i]);
            tmp_fr.mul_assign(aggregation_challenge);
            aggregated_value.add_assign(tmp_fr);
        }
        
        for (uint i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            aggregation_challenge.mul_assign(state.v);

            tmp_fr.assign(proof.permutation_polynomials_at_z[i]);
            tmp_fr.mul_assign(aggregation_challenge);
            aggregated_value.add_assign(tmp_fr);
        }
        
        tmp_fr.assign(proof.grand_product_at_z_omega);
        tmp_fr.mul_assign(state.u);
        aggregated_value.add_assign(tmp_fr);
        
        // final compute
        commitment_aggregation.point_sub_assign(PairingsBn254.P1().point_mul(aggregated_value));
        
        PairingsBn254.G1Point memory pair_with_generator = commitment_aggregation;
        pair_with_generator.point_add_assign(proof.opening_at_z_proof.point_mul(state.z));
        
        tmp_fr.assign(state.z);
        tmp_fr.mul_assign(vk.omega);
        tmp_fr.mul_assign(state.u);
        pair_with_generator.point_add_assign(proof.opening_at_z_omega_proof.point_mul(tmp_fr));
        
        PairingsBn254.G1Point memory pair_with_x = proof.opening_at_z_omega_proof.point_mul(state.u);
        pair_with_x.point_add_assign(proof.opening_at_z_proof);
        pair_with_x.negate();
        
        return PairingsBn254.pairingProd2(pair_with_generator, PairingsBn254.P2(), pair_with_x, vk.g2_x);
    }

    function verify_initial(
        PartialVerifierState memory state, 
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        require(proof.input_values.length == vk.num_inputs);
        require(vk.num_inputs >= 1);
        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.new_transcript();
        
        for (uint256 i = 0; i < proof.wire_commitments.length; i++) {
            transcript.update_with_g1(proof.wire_commitments[i]);
        }
        
        state.beta = transcript.get_challenge(0);
        state.gamma = transcript.get_challenge(1);
        
        transcript.update_with_g1(proof.grand_product_commitment);
        state.alpha = transcript.get_challenge(0);
        
        for (uint256 i = 0; i < proof.quotient_poly_commitments.length; i++) {
            transcript.update_with_g1(proof.quotient_poly_commitments[i]);
        }
    
        state.z = transcript.get_challenge(0);
        
        uint256[] memory lagrange_poly_numbers = new uint256[](vk.num_inputs);
        for (uint256 i = 0; i < lagrange_poly_numbers.length; i++) {
            lagrange_poly_numbers[i] = i;
        }
        
    	state.cached_lagrange_evals = likely_batch_evaluate_lagrange_poly_out_of_domain(
            lagrange_poly_numbers,
            vk.domain_size, 
            vk.omega, state.z
        );

        bool valid = verify_at_z(state, proof, vk);

        if (valid == false) {
            return false;
        }
        
        for (uint256 i = 0; i < proof.wire_values_at_z.length; i++) {
            transcript.update_with_fr(proof.wire_values_at_z[i]);
        }
        
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            transcript.update_with_fr(proof.permutation_polynomials_at_z[i]);
        }
        
        transcript.update_with_fr(proof.grand_product_at_z_omega);
        transcript.update_with_fr(proof.quotient_polynomial_at_z);
        transcript.update_with_fr(proof.linearization_polynomial_at_z);
        
        state.v = transcript.get_challenge(0);
        
        transcript.update_with_g1(proof.opening_at_z_proof);
        transcript.update_with_g1(proof.opening_at_z_omega_proof);
        state.u = transcript.get_challenge(0);

        return true;
    }
    
    // This verifier is for a PLONK with a state width 3
    // and main gate equation
    // q_l(X) * a(X) + 
    // q_r(X) * b(X) + 
    // q_o(X) * c(X) + 
    // q_m(X) * a(X) * b(X) + 
    // q_constants(X)    
    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {
        PartialVerifierState memory state;
        
        bool valid = verify_initial(state, proof, vk);
        
        if (valid == false) {
            return false;
        }
        
        valid = verify_commitments(state, proof, vk);
        
        return valid;
    }
}

contract ConcreteVerifier is Plonk4VerifierWithAccessToDNext { // TODO: maybe we should input all these values by ourselves
    function get_verification_key() internal pure returns(VerificationKey memory vk) {
        vk.domain_size = 4;
        vk.num_inputs = 4;
        vk.omega = PairingsBn254.new_fr(0x30644e72e131a029048b6e193fd841045cea24f6fd736bec231204708f703636);
        vk.selector_commitments[0] = PairingsBn254.new_g1(
            0x2c2c205c4543fbf96d4e00a1867fa9fd6b4adc8e1e63967d6d2bf4080378e16d,
            0x303e034fa6dc771b84fcc4cdaf0e55be886adcf8f8449bfb462e4953bb10e866
        );
        vk.selector_commitments[1] = PairingsBn254.new_g1(
            0x2c2c205c4543fbf96d4e00a1867fa9fd6b4adc8e1e63967d6d2bf4080378e16d,
            0x303e034fa6dc771b84fcc4cdaf0e55be886adcf8f8449bfb462e4953bb10e866 
        );
        vk.selector_commitments[2] = PairingsBn254.new_g1(
            0x1,
            0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd45
        );
        vk.selector_commitments[3] = PairingsBn254.new_g1(
            0xd710d3e1642b7126ca5592fcb6c6e131b9cd45091ed2af3e8ef14d67565e957,
            0x17cb785604229f409460c98c26fe3789dc58569b9dab073a3e6174015e433f3
        );
        vk.selector_commitments[4] = PairingsBn254.new_g1(
            0x0,
            0x0
        );
        
        vk.permutation_commitments[0] = PairingsBn254.new_g1(
            0x1cfd17abd96b518c0552897243def209470f8e12bb32c37e574adbb20c6a900a,
            0xe43639b48ba7e2e8c57482814e90fe6a5e142b17e5ea9aa68a036d2fb0b6a13
        );
        vk.permutation_commitments[1] = PairingsBn254.new_g1(
            0x24f67fdf66438e753c432c5ce23926002c6631f4a18b99ce78c3dee532bcdb71,
            0xc0f9fb65aefe05f3e232a593dd6e68cbda102e2dfc0d8e58501cb5be63befd0 
        );
        vk.permutation_commitments[2] = PairingsBn254.new_g1(
            0x296fcbb175040681d71ad7bc5b3d2a76ba9161c91620b2352527f73947f220f1,
            0x1e89a0a49a0e3262eb9ff0c223bf1b53e5793c4ed5069a80403f0879dcd0a1d
        );
        
        vk.permutation_non_residues[0] = PairingsBn254.new_fr(
            0x0000000000000000000000000000000000000000000000000000000000000005
        );
        vk.permutation_non_residues[1] = PairingsBn254.new_fr(
            0x0000000000000000000000000000000000000000000000000000000000000007
        );
        
        vk.g2_x = PairingsBn254.new_g2(
            [0x19ec5ad0b21f625e7f8eca8841cbd8770106eb7efb06989fd3298ed567597995,
             0x27f29cdb17cbfd46d12d168cadb420d700119e2f04180723025b8be9693e91f6],
            [0x2065616c257a1bfcc4e3ed1175bf0cd77809680ffeb60b274e9d4dad14cbada7,
             0x2423b65cb38edc67d762d937621acc47d1e03744a60ba78c95521ae0064c2463]
        );
    }

	/* public_inputs: [0, 0, 0, 0] */
	/* serialized_proof: 
	["594476778351543693623182262097834998830026499323890699872502437340265883736",
    "11138148742172743821269001527522774419034328252284481412587122967722698849290",
    "16988511013495360245480563426690919669258873728847167534483020190034280476142",
    "7690884383686701489724130432182524583154338824945792131033321852580922834835",
    "14750912893914679006872222158766588469857665657000202622551542322840640583548",
    "7449771727960574255059109935877218648722395938867957119204101734800514363612",
    "1221793087184080053543516341515648316603827204892680018821632993866562284772",
    "1393019551780324211868162829881141593081535791717041811397170102358231997158",
    "106240385429922790006304914823251186024135017389955182129201281796587846931",
    "13622643429880178749639693898343076478307946493431548961689529255788054499889",
    "9192441544880966631260314851905938841390508492840335357157284939279250552077",
    "7570559003449969015057120557838808432517822370579950639299754129070419758843",
    "7087346970662377867758448360849252754250060148089969125635685678506526617607",
    "7395228979775793488657685381967937991224457450131796282771109905859621767169",
    "4541048734516535128845231557055927225042336988981539618395153840322584980389",
    "15776193294210745169799383860962193414359706174313127853847803580994615376591",
    "6584279185101651069614797372532458148942679994682461772209192873084913367960",
    "11319220531976246239681551483899404477295303665041681461820376667619378219900",
    "15794956774158904351606747575091145020757228006229228194018827461355852755859",
    "10226390296738572000586020641836306372140448879958659378103364447108688106024",
    "18292239692626280844447445246203069344773166963401944699579184196610619976184",
    "2545119322393522288860947374040065139955988451279616812980961106510828092291",
    "18889419035524852273906886002593758661088331730692952604750294729011358113690",
    "7211031983153934433068056236562777502906276950910502670498850085188670954815",
    "8723629029486532637554723971251827809244379053483711686717418536067147115884",
    "8101898248330553924390358101659418761457918957779251740823643841466732107162"]
	*/

    function deserialize_proof(
        uint256 expected_inputs, 
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) internal pure returns(Proof memory proof) {
        assert(expected_inputs == public_inputs.length);
        proof.input_values = new uint256[](expected_inputs);
        for (uint256 i = 0; i < expected_inputs; i++) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j = 0;
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_commitments[i] = PairingsBn254.new_g1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.grand_product_commitment = PairingsBn254.new_g1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.quotient_poly_commitments[i] = PairingsBn254.new_g1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_values_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.grand_product_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.quotient_polynomial_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;

        proof.linearization_polynomial_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;
    
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            proof.permutation_polynomials_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }

        proof.opening_at_z_proof = PairingsBn254.new_g1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;

        proof.opening_at_z_omega_proof = PairingsBn254.new_g1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
    
        j += 2;
        assert(j == serialized_proof.length);
    }
    
    function verify(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        uint256 expected_inputs = vk.num_inputs;

        Proof memory proof = deserialize_proof(expected_inputs, public_inputs, serialized_proof);

        bool valid = verify(proof, vk);

        return valid;
    }  
}
