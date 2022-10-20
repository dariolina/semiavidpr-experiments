use super::*;


use ark_bls12_381::{Bls12_381};
use ark_bn254::{Bn254};


#[test]
fn test_kzg_commit_bls12_381() {
    _test_kzg_commit::<Bls12_381>()
}

#[test]
fn test_kzg_commit_bn254() {
    _test_kzg_commit::<Bn254>()
}

fn _test_kzg_commit<E: PairingEngine>() {
    let mut rng = ark_std::rand::thread_rng();
    
    let scheme = SemiAvidPr::<E>::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    for i in 0..scheme.k {
        // internal method
        let commitment1 = scheme.commit_column(&data_uncoded, i);

        // explicit method
        let poly_evals = Evaluations::from_vec_and_domain(data_uncoded.iter().map(|r| r[i]).collect(), scheme.domain_polycommit);
        let poly_poly = poly_evals.interpolate();
        let commitment2 = SemiAvidPr::unwrap_commitment(KZG10::commit(&scheme.kzg10_ck, &poly_poly, None, None).unwrap());
        
        // compute KZG commitment manually to double-check calculation
        let mut commitment3 = E::G1Projective::zero();
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
    _test_kzg_commitment_unwrap_wrap::<Bls12_381>()
}

#[test]
fn test_kzg_commitment_unwrap_wrap_bn254() {
    _test_kzg_commitment_unwrap_wrap::<Bn254>()
}

fn _test_kzg_commitment_unwrap_wrap<E: PairingEngine>() {
    let mut rng = ark_std::rand::thread_rng();
    
    let scheme = SemiAvidPr::<E>::setup(&mut rng, 2, 1, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);

    let poly_evals = Evaluations::from_vec_and_domain(data_uncoded.iter().map(|r| r[0]).collect(), scheme.domain_polycommit);
    let poly_poly = poly_evals.interpolate();
    let commitment = KZG10::commit(&scheme.kzg10_ck, &poly_poly, None, None).unwrap();

    assert_eq!(commitment, SemiAvidPr::wrap_commitment(SemiAvidPr::unwrap_commitment(commitment.clone())));
}

#[test]
fn test_filesizes() {
    let mut rng = ark_std::rand::thread_rng();
    let scheme = SemiAvidPr::<Bls12_381>::setup(&mut rng, 512, 256, 1024);
    assert_eq!(scheme.get_filesize_in_bytes(), 254 * 256*1024 / 8);
    let scheme = SemiAvidPr::<Bn254>::setup(&mut rng, 512, 256, 1024);
    assert_eq!(scheme.get_filesize_in_bytes(), 253 * 256*1024 / 8);
}
#[test]
fn test_commit_commit_bls12_381() {
    _test_commit_commit::<Bls12_381>()
}
//tests 'commitment to commitments' i.e a polynomial where commitments are coefficients evaluated at public parameters
fn _test_commit_commit<E: PairingEngine>() {
    let mut rng = ark_std::rand::thread_rng();
    
    let scheme = SemiAvidPr::<E>::setup(&mut rng, 16, 8, 1024);
    let data_uncoded = scheme.generate_random_file(&mut rng);
    let mut comm_root = E::G1Projective::zero();

    //i iterates columns
    for i in 0..scheme.k {
        let poly_evals = Evaluations::from_vec_and_domain(data_uncoded.iter().map(|r| r[i]).collect(), scheme.domain_polycommit);
        let poly_poly = poly_evals.interpolate();

        //compute corresponding term into commitment root
        let mut comm_root_term = E::G1Projective::zero();

        assert_eq!(scheme.kzg10_ck.powers_of_g.len(), poly_poly.coeffs.len()+scheme.k);

        //j iterates rows
        for j in 0..scheme.kzg10_ck.powers_of_g.len()-scheme.k {
            comm_root_term += scheme.kzg10_ck.powers_of_g[j+i].mul(poly_poly.coeffs[j]);
        }
        
        comm_root += comm_root_term;
        
    }   

    //should be anything but 0   
    assert!(
            !comm_root.is_zero(),
            "oops"
     );
}