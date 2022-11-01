use crate::SemiAvidPr;
use ark_bls12_381::G1Projective;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::Zero;
use ark_poly::Evaluations;
use ark_poly_commit::kzg10::KZG10;

#[test]
fn test_kzg_commit_bls12_381() {
    let mut rng = rand::thread_rng();

    let scheme = SemiAvidPr::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    for i in 0..scheme.k {
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
    for i in 0..scheme.k {
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
            poly_poly.coeffs.len() + scheme.k
        );

        //j iterates rows
        for j in 0..scheme.kzg10_ck.powers_of_g.len() - scheme.k {
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
    assert_eq!(data_coded[0][scheme.k - 1], data_uncoded[0][scheme.k - 1]);

    assert!(scheme.disperse_verify_chunks_systematic(&source_column_commitments, &data_coded));
}
