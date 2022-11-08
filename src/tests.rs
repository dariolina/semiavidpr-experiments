use crate::SemiAvidPr;
use ark_bls12_381::{G1Projective, Fr};
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{Zero, FftField};
use ark_poly::{evaluations::univariate::Evaluations, 
    polynomial::univariate::DensePolynomial,Polynomial,UVPolynomial,
    EvaluationDomain,Radix2EvaluationDomain};
use ark_poly_commit::kzg10::KZG10;

#[test]
fn test_kzg_commit_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    for i in 0..scheme.uncoded_chunks {
        // internal method
        let commitment1 = scheme.commit_column(&data_uncoded, i);

        // explicit method
        let poly_evals = Evaluations::from_vec_and_domain(
            data_uncoded.iter().map(|r| r[i]).collect(),
            scheme.domain_polycommit,
        );
        let poly_poly = poly_evals.interpolate();
        let commitment2 = SemiAvidPr::unwrap_commitment(
            KZG10::commit(&scheme.kzg10_ck, &poly_poly, None, None).unwrap(),
        );

        // compute KZG commitment manually to double-check calculation
        let mut commitment3 = G1Projective::zero();
        assert_eq!(scheme.kzg10_ck.powers_of_g.len(), poly_poly.coeffs.len());
        for j in 0..scheme.kzg10_ck.powers_of_g.len() {
            commitment3 += scheme.kzg10_ck.powers_of_g[j].mul(poly_poly.coeffs[j]);
        }
        let commitment3 = commitment3.into_affine();

        // the above three better be equal
        assert_eq!(commitment1, commitment2);
        assert_eq!(commitment1, commitment3);
        assert_eq!(commitment2, commitment3);
    }
}

#[test]
fn test_kzg_commitment_unwrap_wrap_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 2, 1, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    let poly_evals = Evaluations::from_vec_and_domain(
        data_uncoded.iter().map(|r| r[0]).collect(),
        scheme.domain_polycommit,
    );
    let poly_poly = poly_evals.interpolate();
    let commitment = KZG10::commit(&scheme.kzg10_ck, &poly_poly, None, None).unwrap();

    assert_eq!(
        commitment,
        SemiAvidPr::wrap_commitment(SemiAvidPr::unwrap_commitment(commitment.clone()))
    );
}

#[test]
fn test_filesizes() {
    let mut rng = rand::thread_rng();
    let scheme = SemiAvidPr::setup(&mut rng, 512, 256, 1024);
    assert_eq!(scheme.get_filesize_in_bytes(), 254 * 256 * 1024 / 8);
}

//tests 'commitment to commitments' i.e a polynomial where commitments are coefficients evaluated at public parameters
#[test]
fn test_commit_commit1_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);
    let mut comms_list = Vec::new();
    let mut comm_root1 = G1Projective::zero();

    let mut comm_root2 = G1Projective::zero();

    //i iterates columns
    for i in 0..scheme.uncoded_chunks {
        let poly_evals = Evaluations::from_vec_and_domain(
            data_uncoded.iter().map(|r| r[i]).collect(),
            scheme.domain_polycommit,
        );
        let poly_poly = poly_evals.interpolate();

        //compute commitment
        let mut commitment = G1Projective::zero();
        //compute corresponding term into commitment root
        let mut comm_root_term1 = G1Projective::zero();
        let mut comm_root_term2 = G1Projective::zero();

        assert_eq!(
            scheme.kzg10_ck.powers_of_g.len(),
            poly_poly.coeffs.len() + scheme.uncoded_chunks
        );

        //j iterates rows
        for j in 0..scheme.kzg10_ck.powers_of_g.len() - scheme.uncoded_chunks {
            commitment += scheme.kzg10_ck.powers_of_g[j].mul(poly_poly.coeffs[j]);
            comm_root_term1 += scheme.kzg10_ck.powers_of_g[j + i].mul(poly_poly.coeffs[j]);
            comm_root_term2 += commitment + scheme.kzg10_ck.powers_of_g[i].mul(poly_poly.coeffs[j]);
        }
        comms_list.push(commitment);
        comm_root1 += comm_root_term1;
        comm_root2 += comm_root_term2;
    }
    assert_ne!(comm_root1, comm_root2);
}

//tests systematic encoding of data and commitments
#[test]
fn test_systematic_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    let source_column_commitments = scheme.disperse_compute_column_commitments(&data_uncoded);

    let data_coded = scheme.disperse_encode_rows_lagrange(&data_uncoded);
    //assert uncoded data matches the first indices of coded
    assert_eq!(data_coded[0][0], data_uncoded[0][0]);
    assert_eq!(
        data_coded[0][scheme.uncoded_chunks - 1],
        data_uncoded[0][scheme.uncoded_chunks - 1]
    );

    //assert coded data doesn't copy uncoded
    assert_ne!(data_coded[0][scheme.uncoded_chunks], data_coded[0][0]);

    assert!(scheme.disperse_verify_chunks_systematic(&source_column_commitments, &data_coded));
}
#[test]
fn test_fft_domains() {
    // Tests that the ffts output the correct result.
    // This assumes a correct polynomial evaluation at point procedure.
    // It tests consistency of FFT/IFFT, and coset_fft/coset_ifft,
    // along with testing that each individual evaluation is correct.

    // Runs in time O(degree^2)
    let mut rng = rand::thread_rng();

    let uncoded_chunks = 8;
    let coded_chunks = 16;

    let scheme = SemiAvidPr::setup(&mut rng, coded_chunks, uncoded_chunks, 32);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    //make some random data
    let domain = Radix2EvaluationDomain::<Fr>::new(uncoded_chunks).unwrap();
         
    let source_evals = Evaluations::from_vec_and_domain(
                data_uncoded[0].iter().copied().collect(),
                domain,
            );

    let poly_poly = source_evals.interpolate();
    
    let poly_coset_evals = domain.coset_fft(&poly_poly.coeffs);
    
    for (i, x) in domain.elements().enumerate() {
            let coset_x = Fr::multiplicative_generator() * x;

            assert_eq!(data_uncoded[0][i], poly_poly.evaluate(&x));
            assert_eq!(poly_coset_evals[i], poly_poly.evaluate(&coset_x));
            //coset evals are not the same as source
            assert_ne!(data_uncoded[0][i], poly_coset_evals[i]);
        }

        let poly_from_coset =
            DensePolynomial::from_coefficients_vec(domain.coset_ifft(&poly_coset_evals));

        assert_eq!(
            poly_poly, poly_from_coset,
            "degree = {}, domain size = {}",
            poly_poly.degree(), uncoded_chunks
        );
    
}
#[test]
fn test_fft_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    let source_column_commitments = scheme.disperse_compute_column_commitments(&data_uncoded);

    let data_coded = scheme.disperse_encode_rows_fft(&data_uncoded);
    //assert uncoded data matches the even indices of coded
    assert_eq!(data_coded[0][0], data_uncoded[0][0]);
    assert_eq!(data_coded[0][1], data_uncoded[0][1]);

    //assert coded data doesn't copy uncoded
    assert_ne!(data_coded[0][scheme.uncoded_chunks], data_coded[0][0]);

    assert!(scheme.disperse_verify_chunks_fft(&source_column_commitments, &data_coded));
}

#[test]
fn prototype_encoding_decoding() {
    let mut rng = rand::thread_rng();

    let uncoded_chunks = 8;
    let coded_chunks = 16;

    let scheme = SemiAvidPr::setup(&mut rng, coded_chunks, uncoded_chunks, 32);

    let data_uncoded = scheme.generate_random_file(&mut rng);

    let data_coded = scheme.disperse_encode_rows_fft(&data_uncoded);
    // Chunks from which erasure decoding will happen, take last `uncoded_chunks`
    let coded_entries: Vec<usize> = (0..coded_chunks).rev().take(uncoded_chunks).rev().collect();

    let data_coded_downloaded = scheme.retrieve_download_chunks(&data_coded, &coded_entries);

    let data_uncoded_downloaded = scheme.retrieve_decode_rows_fft(&data_coded_downloaded);

    assert_eq!(data_uncoded, data_uncoded_downloaded);
}

#[test]
fn prototype_commit_prove_verify() {
    let mut rng = rand::thread_rng();

    let uncoded_chunks = 8;
    let coded_chunks = 16;

    let scheme = SemiAvidPr::setup(&mut rng, coded_chunks, uncoded_chunks, 32);

    let data = scheme
        .generate_random_file(&mut rng)
        .into_iter()
        .next()
        .unwrap();

    let polynomial = scheme.create_polynomial(data.clone());
    let commitment = scheme.commit(&polynomial);

    for (idx, d) in data.into_iter().enumerate() {
        let witness_commitment = scheme.commit(&scheme.create_witness_polynomial(&polynomial, idx));

        assert!(
            scheme.verify(d, idx, commitment, witness_commitment),
            "Idx {idx}"
        );
    }
}
