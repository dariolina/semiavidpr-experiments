use ark_ff::fields::{Field, PrimeField, FpParameters};
use ark_ec::{ProjectiveCurve, AffineCurve, PairingEngine};
use ark_poly::{
    Polynomial, UVPolynomial,
    EvaluationDomain, GeneralEvaluationDomain,
    polynomial::univariate::{DensePolynomial},
    evaluations::univariate::{Evaluations},
};
use ark_poly_commit::{
    PCRandomness,
    kzg10::{KZG10, Powers, VerifierKey, Commitment, Proof, Randomness},
};
use ark_std::{Zero, One, UniformRand, start_timer, end_timer};

use rand::{Rng};


mod utils;
use crate::utils::{Matrix};

#[cfg(test)]
mod tests;


#[allow(non_snake_case)]
pub struct SemiAvidPr<'a, E: PairingEngine> {
    n: usize,
    k: usize,
    L: usize,

    domain_polycommit: GeneralEvaluationDomain<E::Fr>,
    domain_encoding: GeneralEvaluationDomain<E::Fr>,

    kzg10_ck: Powers<'a, E>,
    kzg10_vk: VerifierKey<E>,
}


impl<E: PairingEngine> SemiAvidPr<'_, E> {
    #[allow(non_snake_case)]
    pub fn setup<R: Rng + ?Sized>(mut rng: &mut R, n: usize, k: usize, L: usize) -> Self {
        assert!(n.is_power_of_two());
        assert!(L.is_power_of_two());

        let timer = start_timer!(|| "Creating evaluation domains");
        let domain_polycommit: GeneralEvaluationDomain<E::Fr> = ark_poly::domain::EvaluationDomain::<E::Fr>::new(L).unwrap();
        let domain_encoding: GeneralEvaluationDomain<E::Fr> = ark_poly::domain::EvaluationDomain::<E::Fr>::new(n).unwrap();
        end_timer!(timer);

        let timer = start_timer!(|| "KZG setup and preprocessing of setup");
        let kzg10_pp = KZG10::<E, DensePolynomial<E::Fr>>::setup(L+k-1, false, &mut rng).unwrap();

        // https://github.com/arkworks-rs/poly-commit/blob/4d78d534cb55a9b13f34dd76b9702cae3ab2a2a1/src/kzg10/mod.rs#L459
        let (kzg10_ck, kzg10_vk) = {
            let powers_of_g = kzg10_pp.powers_of_g[..=(L+k-1)].to_vec();
            let powers_of_gamma_g = (0..=(L+k-1))
                .map(|i| kzg10_pp.powers_of_gamma_g[&i])
                .collect();

            let powers = Powers {
                powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
                powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
            };
            let vk = VerifierKey::<E> {
                g: kzg10_pp.powers_of_g[0],
                gamma_g: kzg10_pp.powers_of_gamma_g[&0],
                h: kzg10_pp.h,
                beta_h: kzg10_pp.beta_h,
                prepared_h: kzg10_pp.prepared_h.clone(),
                prepared_beta_h: kzg10_pp.prepared_beta_h.clone(),
            };

            (powers, vk)
        };
        end_timer!(timer);

        Self {
            n, k, L,

            domain_polycommit,
            domain_encoding,

            kzg10_ck,
            kzg10_vk,
        }
    }


    pub fn get_filesize(&self) -> usize {
        (<E::Fr as PrimeField>::Params::CAPACITY as usize) * self.k * self.L
    }

    pub fn get_filesize_in_bytes(&self) -> u64 {
        (self.get_filesize() / 8) as u64
    }

    pub fn get_num_column_commitments(&self) -> usize {
        return self.k
    }

    pub fn get_num_row_encodings(&self) -> usize {
        return self.L
    }

    pub fn get_num_chunk_verifications(&self) -> usize {
        return self.n
    }

    pub fn get_num_downloaded_chunk_verifications(&self) -> usize {
        return self.k
    }

    pub fn get_num_row_decodings(&self) -> usize {
        return self.L
    }


    pub fn generate_random_file<R: Rng + ?Sized>(&self, mut rng: &mut R) -> Vec<Vec<E::Fr>> {
        let mut data = vec![vec![E::Fr::zero(); self.k]; self.L];

        let timer = start_timer!(|| "Sampling random field elements");
        for i in 0..self.k {
            for j in 0..self.L {
                data[j][i] = E::Fr::rand(&mut rng);
            }
        }
        end_timer!(timer);

        data
    }


    fn unwrap_commitment(c: (Commitment<E>, Randomness<E::Fr, DensePolynomial<E::Fr>>)) -> E::G1Affine {
        c.0.0
    }

    fn wrap_commitment(c: E::G1Affine) -> (Commitment<E>, Randomness<E::Fr, DensePolynomial<E::Fr>>) {
        (Commitment::<E>(c), Randomness::<E::Fr, DensePolynomial<E::Fr>>::empty())
    }


    fn commit_column(&self, data: &Vec<Vec<E::Fr>>, idx: usize) -> E::G1Affine {
        let timer = start_timer!(|| "Poly evaluations and interpolation");
        let poly_evals = Evaluations::from_vec_and_domain(data.iter().map(|r| r[idx]).collect(), self.domain_polycommit);
        let poly_poly = poly_evals.interpolate();
        end_timer!(timer);
        
        let timer = start_timer!(|| "KZG commitment");
        let commitment = Self::unwrap_commitment(KZG10::commit(&self.kzg10_ck, &poly_poly, None, None).unwrap());
        end_timer!(timer);

        commitment
    }


    fn encode_commitments(&self, column_commitments: &Vec<E::G1Affine>, idx: usize) -> E::G1Affine {
        let timer = start_timer!(|| "'Encoding' of KZG column commitments");
        let mut commitment = E::G1Projective::zero();
        for j in 0..self.k {
            let j_in_field = E::Fr::from_le_bytes_mod_order(&j.to_le_bytes());
            let eval_exponent = self.domain_encoding.element(idx).pow(j_in_field.into_repr());
            commitment += column_commitments[j].mul(eval_exponent);
        }
        let commitment = commitment.into_affine();
        end_timer!(timer);

        commitment
    }
    fn lagrange(&self, i:usize, idx:usize)->E::Fr{
        let mut coef =E::Fr::one();

        let domain_uncoded: GeneralEvaluationDomain<E::Fr> = ark_poly::domain::EvaluationDomain::<E::Fr>::new(self.k).unwrap();
        //assert first element is the field 1
        assert_eq!(domain_uncoded.element(0), E::Fr::one(),"one is not one");

        for j in 0..self.k{
            if j==i{
                continue;
            }
            coef = coef * (domain_uncoded.element(idx)-domain_uncoded.element(j))/(domain_uncoded.element(i)-domain_uncoded.element(j));
        }
        if idx==i{
            //assert definition of lagrange basis poly
            assert_eq!(coef, E::Fr::one());
        }
        coef
    }

    fn encode_commitments_systematic(&self, column_commitments: &Vec<E::G1Affine>, idx: usize) -> E::G1Affine {
        let timer = start_timer!(|| "'Encoding' of KZG column commitments");
        let mut commitment = E::G1Projective::zero();
        for j in 0..self.k {
            let coef = self.lagrange(j, idx);
            commitment += column_commitments[j].mul(coef);
        }
        let commitment = commitment.into_affine();
        end_timer!(timer);

        commitment
    }
    
    

    pub fn disperse_compute_column_commitments(&self, data_uncoded: &Vec<Vec<E::Fr>>) -> Vec<E::G1Affine> {
        let mut column_commitments = Vec::new();

        let timer_outer = start_timer!(|| "Computing column commitments");
        for i in 0..self.k {
            let timer_inner = start_timer!(|| format!("Column {}", i));

            let commitment = self.commit_column(&data_uncoded, i);
            column_commitments.push(commitment);

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        column_commitments
    }


    pub fn disperse_encode_rows(&self, data_uncoded: &Vec<Vec<E::Fr>>) -> Vec<Vec<E::Fr>> {
        let mut data_coded = Vec::new();

        let timer_outer = start_timer!(|| "Encoding rows");
        for j in 0..self.L {
            let timer_inner = start_timer!(|| format!("Row {}", j));

            let poly_poly = DensePolynomial::<E::Fr>::from_coefficients_slice(&data_uncoded[j]);
            let poly_evals = poly_poly.evaluate_over_domain(self.domain_encoding);
            data_coded.push(poly_evals.evals);

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        data_coded
    }

    pub fn disperse_encode_rows_systematic(&self, data_uncoded: &Vec<Vec<E::Fr>>) -> Vec<Vec<E::Fr>> {
        let mut data_coded = Vec::new();

        let timer_outer = start_timer!(|| "Encoding rows");
        let domain_uncoded: GeneralEvaluationDomain<E::Fr> = ark_poly::domain::EvaluationDomain::<E::Fr>::new(self.k).unwrap();
        for j in 0..self.L {
            let timer_inner = start_timer!(|| format!("Row {}", j));
            //assert expected length of source data
            assert_eq!(data_uncoded[j].len(), self.k);

            let mut poly_evals = Evaluations::from_vec_and_domain(data_uncoded[j].iter().copied().collect(), domain_uncoded);
            let poly_poly = poly_evals.interpolate();
            //assert expected degree of interpolated polynomial
            assert_eq!(poly_poly.degree(), self.k-1);

            poly_evals = poly_poly.evaluate_over_domain(self.domain_encoding);

            data_coded.push(poly_evals.evals);
            assert_eq!(data_coded[j].len(), self.n);
            

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        data_coded
    }

    pub fn disperse_verify_chunks(&self, column_commitments: &Vec<E::G1Affine>, data_coded: &Vec<Vec<E::Fr>>) -> bool {
        let timer_outer = start_timer!(|| "Checking coded columns");
        for i in 0..self.n {
            let timer_inner = start_timer!(|| format!("Column {}", i));

            let commitment = self.commit_column(&data_coded, i);
            let commitment_check = self.encode_commitments(&column_commitments, i);
            

            if commitment != commitment_check {
                return false;
            }

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        true
    }
    pub fn disperse_verify_chunks_systematic(&self, column_commitments: &Vec<E::G1Affine>, data_coded: &Vec<Vec<E::Fr>>) -> bool {
        let timer_outer = start_timer!(|| "Checking coded columns");
        for i in 0..self.n {
            let timer_inner = start_timer!(|| format!("Column {}", i));

            let commitment = self.commit_column(&data_coded, i);
            let commitment_check = self.encode_commitments_systematic(&column_commitments, i);
            if commitment != commitment_check {
                return false;
            }

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        true
    }

    pub fn retrieve_download_chunks(&self, data_coded: &Vec<Vec<E::Fr>>, idxs_download_nodes: &Vec<usize>) -> Vec<Vec<E::Fr>> {
        let timer = start_timer!(|| "Downloading chunks");
        let data_coded_downloaded = data_coded.iter().map(|r| idxs_download_nodes.iter().map(|&i| r[i]).collect()).collect();
        end_timer!(timer);

        data_coded_downloaded
    }


    pub fn retrieve_verify_chunks(&self, column_commitments: &Vec<E::G1Affine>, data_coded_downloaded: &Vec<Vec<E::Fr>>, idxs_download_nodes: &Vec<usize>) -> bool {
        let timer_all = start_timer!(|| "Verifying downloaded chunks");

        let timer_encode_commitments = start_timer!(|| "'Encoding' of column commitments to coded chunk commitments");
        let mut column_commitments_projective: Vec<E::G1Projective> = column_commitments.iter().map(|h| h.clone().into()).collect();
        column_commitments_projective.resize(self.domain_encoding.size(), E::G1Projective::zero());
        self.domain_encoding.fft_in_place(&mut column_commitments_projective);
        let coded_chunk_commitments_affine: Vec<E::G1Affine> = column_commitments_projective.iter().map(|h| h.into_affine()).collect();
        end_timer!(timer_encode_commitments);

        let timer_outer = start_timer!(|| "Checking downloaded coded columns");
        for (idx, col) in idxs_download_nodes.iter().enumerate() {
            let timer_inner = start_timer!(|| format!("Column {}", idx));

            let commitment = self.commit_column(&data_coded_downloaded, idx);

            // let commitment_check = self.encode_commitments(&column_commitments, *col);
            // if commitment != commitment_check {
            //     return false;
            // }

            if commitment != coded_chunk_commitments_affine[*col] {
                return false;
            }

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        end_timer!(timer_all);
        true
    }
    

    pub fn retrieve_prepare_decoding(&self, idxs_download_nodes: &Vec<usize>) -> Matrix<E::Fr> {
        assert!(idxs_download_nodes.len() == self.k);

        let mut matrix = Vec::new();
        for i in 0..self.k {
            let i_in_field = E::Fr::from_le_bytes_mod_order(&i.to_le_bytes());
            matrix.push(idxs_download_nodes.iter().map(|&j| self.domain_encoding.element(j).pow(i_in_field.into_repr())).collect());
        }

        Matrix::from_nested_vec(self.k, self.k, matrix).invert()
    }

    
    pub fn retrieve_decode_rows(&self, data_coded_downloaded: &Vec<Vec<E::Fr>>, decoder_aux: &Matrix<E::Fr>) -> Vec<Vec<E::Fr>> {
        assert!(decoder_aux.height() == decoder_aux.width());
        let mut data_decoded = Vec::new();

        let timer_outer = start_timer!(|| "Decoding rows");
        for j in 0..self.L {
            let timer_inner = start_timer!(|| format!("Row {}", j));

            assert!(data_coded_downloaded[j].len() == decoder_aux.height());
            data_decoded.push((0..decoder_aux.width()).map(|col| (0..decoder_aux.height()).map(|row| decoder_aux.get(row, col) * data_coded_downloaded[j][row]).sum::<E::Fr>()).collect());

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        data_decoded
    }


    pub fn sampling_open_entry(&self, column_commitments: &Vec<E::G1Affine>, data_uncoded: &Vec<Vec<E::Fr>>, row: usize, col: usize) -> (E::Fr, usize, usize, Vec<E::G1Affine>, Proof<E>) {
        let timer = start_timer!(|| "Poly evaluations and interpolation");
        let poly_evals = Evaluations::from_vec_and_domain(data_uncoded.iter().map(|r| r[col]).collect(), self.domain_polycommit);
        let poly_poly = poly_evals.interpolate();
        end_timer!(timer);

        let timer = start_timer!(|| "KZG proof");
        // let proof = KZG10::open(&self.kzg10_ck, &poly_poly, self.domain_polycommit.element(row), None).unwrap();
        // Unfortunately, KZG10::open() is pub(crate) only, so inline ... >>>
        let point = self.domain_polycommit.element(row);
        assert!(poly_poly.degree() + 1 <= self.kzg10_ck.size());
        let divisor = DensePolynomial::<E::Fr>::from_coefficients_vec(vec![-point, E::Fr::one()]);
        let witness_polynomial = &poly_poly / &divisor;
        assert!(witness_polynomial.degree() + 1 <= self.kzg10_ck.size());
        let proof = Self::unwrap_commitment(KZG10::commit(&self.kzg10_ck, &witness_polynomial, None, None).unwrap());
        let proof = Proof { w: proof, random_v: None };
        // <<< ... end of inline!
        end_timer!(timer);

        (data_uncoded[row][col], row, col, column_commitments.clone(), proof)
    }


    pub fn sampling_verify_entry(&self, (value, row, col, column_commitments, proof): (E::Fr, usize, usize, Vec<E::G1Affine>, Proof<E>)) -> bool {
        let timer = start_timer!(|| "KZG check");
        let commitment = Self::wrap_commitment(column_commitments[col]).0;
        let point = self.domain_polycommit.element(row);
        let ret_val = KZG10::<E, DensePolynomial::<E::Fr>>::check(&self.kzg10_vk, &commitment, point, value, &proof).unwrap();
        end_timer!(timer);
        
        ret_val
    }
}
