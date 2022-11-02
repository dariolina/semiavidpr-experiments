use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective};
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::fields::{Field, FpParameters, PrimeField};
use ark_poly::{
    evaluations::univariate::Evaluations, polynomial::univariate::DensePolynomial,
    EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_poly_commit::{
    kzg10::{Commitment, Powers, Proof, Randomness, VerifierKey, KZG10},
    PCRandomness,
};
use ark_std::{end_timer, start_timer, One, UniformRand, Zero};

use rand::Rng;

mod utils;
use crate::utils::Matrix;

#[cfg(test)]
mod tests;

pub struct SemiAvidPr<'a> {
    coded_chunks: usize,
    uncoded_chunks: usize,
    chunk_length: usize,

    domain_polycommit: GeneralEvaluationDomain<Fr>,
    domain_encoding: GeneralEvaluationDomain<Fr>,

    kzg10_ck: Powers<'a, Bls12_381>,
    kzg10_vk: VerifierKey<Bls12_381>,
}

impl SemiAvidPr<'_> {
    pub fn setup<R: Rng + ?Sized>(
        mut rng: &mut R,
        coded_chunks: usize,
        uncoded_chunks: usize,
        chunk_length: usize,
    ) -> Self {
        assert!(coded_chunks.is_power_of_two());
        assert!(chunk_length.is_power_of_two());

        let timer = start_timer!(|| "Creating evaluation domains");
        let domain_polycommit: GeneralEvaluationDomain<Fr> =
            ark_poly::domain::EvaluationDomain::<Fr>::new(chunk_length).unwrap();
        let domain_encoding: GeneralEvaluationDomain<Fr> =
            ark_poly::domain::EvaluationDomain::<Fr>::new(coded_chunks).unwrap();
        end_timer!(timer);

        let timer = start_timer!(|| "KZG setup and preprocessing of setup");
        let kzg10_pp = KZG10::<Bls12_381, DensePolynomial<Fr>>::setup(
            chunk_length + uncoded_chunks - 1,
            false,
            &mut rng,
        )
        .unwrap();

        // https://github.com/arkworks-rs/poly-commit/blob/4d78d534cb55a9b13f34dd76b9702cae3ab2a2a1/src/kzg10/mod.rs#L459
        let (kzg10_ck, kzg10_vk) = {
            let powers_of_g = kzg10_pp.powers_of_g[..=(chunk_length + uncoded_chunks - 1)].to_vec();
            let powers_of_gamma_g = (0..=(chunk_length + uncoded_chunks - 1))
                .map(|i| kzg10_pp.powers_of_gamma_g[&i])
                .collect();

            let powers = Powers {
                powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
                powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
            };
            let vk = VerifierKey::<Bls12_381> {
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
            coded_chunks,
            uncoded_chunks,
            chunk_length,

            domain_polycommit,
            domain_encoding,

            kzg10_ck,
            kzg10_vk,
        }
    }

    pub fn get_filesize(&self) -> usize {
        (<Fr as PrimeField>::Params::CAPACITY as usize) * self.uncoded_chunks * self.chunk_length
    }

    pub fn get_filesize_in_bytes(&self) -> u64 {
        (self.get_filesize() / 8) as u64
    }

    pub fn get_num_column_commitments(&self) -> usize {
        return self.uncoded_chunks;
    }

    pub fn get_num_row_encodings(&self) -> usize {
        return self.chunk_length;
    }

    pub fn get_num_chunk_verifications(&self) -> usize {
        return self.coded_chunks;
    }

    pub fn get_num_downloaded_chunk_verifications(&self) -> usize {
        return self.uncoded_chunks;
    }

    pub fn get_num_row_decodings(&self) -> usize {
        return self.chunk_length;
    }

    pub fn generate_random_file<R: Rng + ?Sized>(&self, mut rng: &mut R) -> Vec<Vec<Fr>> {
        let mut data = vec![vec![Fr::zero(); self.uncoded_chunks]; self.chunk_length];

        let timer = start_timer!(|| "Sampling random field elements");
        for i in 0..self.uncoded_chunks {
            for j in 0..self.chunk_length {
                data[j][i] = Fr::rand(&mut rng);
            }
        }
        end_timer!(timer);

        data
    }

    fn unwrap_commitment(
        c: (Commitment<Bls12_381>, Randomness<Fr, DensePolynomial<Fr>>),
    ) -> G1Affine {
        c.0 .0
    }

    fn wrap_commitment(
        c: G1Affine,
    ) -> (Commitment<Bls12_381>, Randomness<Fr, DensePolynomial<Fr>>) {
        (
            Commitment::<Bls12_381>(c),
            Randomness::<Fr, DensePolynomial<Fr>>::empty(),
        )
    }

    fn commit_column(&self, data: &Vec<Vec<Fr>>, idx: usize) -> G1Affine {
        let timer = start_timer!(|| "Poly evaluations and interpolation");
        let poly_evals = Evaluations::from_vec_and_domain(
            data.iter().map(|r| r[idx]).collect(),
            self.domain_polycommit,
        );
        let poly_poly = poly_evals.interpolate();
        end_timer!(timer);

        let timer = start_timer!(|| "KZG commitment");
        let commitment =
            Self::unwrap_commitment(KZG10::commit(&self.kzg10_ck, &poly_poly, None, None).unwrap());
        end_timer!(timer);

        commitment
    }

    // fn encode_commitments(&self, column_commitments: &Vec<G1Affine>, idx: usize) -> G1Affine {
    //     let timer = start_timer!(|| "'Encoding' of KZG column commitments");
    //     let mut commitment = G1Projective::zero();
    //     for j in 0..self.k {
    //         let j_in_field = Fr::from_le_bytes_mod_order(&j.to_le_bytes());
    //         let eval_exponent = self
    //             .domain_encoding
    //             .element(idx)
    //             .pow(j_in_field.into_repr());
    //         commitment += column_commitments[j].mul(eval_exponent);
    //     }
    //     let commitment = commitment.into_affine();
    //     end_timer!(timer);
    //
    //     commitment
    // }

    fn evaluate_barycentric(&self, j: usize, idx: usize) -> Fr {
        let d = self.uncoded_chunks - 1;
        let d_in_field = Fr::from_le_bytes_mod_order(&d.to_le_bytes());
        let bary_coef = (self
            .domain_encoding
            .element(idx)
            .pow(d_in_field.into_repr())
            - self.domain_encoding.element(1))
            / self.domain_encoding.element(d);

        (self.domain_encoding.element(j)
            / (self.domain_encoding.element(idx) - self.domain_encoding.element(j)))
            * bary_coef
    }

    fn encode_commitments_systematic(
        &self,
        column_commitments: &Vec<G1Affine>,
        idx: usize,
    ) -> G1Affine {
        let timer = start_timer!(|| "'Encoding' of KZG column commitments");
        let mut commitment = G1Projective::zero();

        if idx >= self.uncoded_chunks {
            for j in 0..self.uncoded_chunks {
                let coef = self.evaluate_barycentric(j, idx);
                commitment += column_commitments[j].mul(coef);
            }
        } else {
            commitment += column_commitments[idx].mul(Fr::one());
        }

        end_timer!(timer);
        let commitment = commitment.into_affine();
        commitment
    }

    pub fn disperse_compute_column_commitments(
        &self,
        data_uncoded: &Vec<Vec<Fr>>,
    ) -> Vec<G1Affine> {
        let mut column_commitments = Vec::new();

        let timer_outer = start_timer!(|| "Computing column commitments");
        for i in 0..self.uncoded_chunks {
            let timer_inner = start_timer!(|| format!("Column {}", i));

            let commitment = self.commit_column(&data_uncoded, i);
            column_commitments.push(commitment);

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        column_commitments
    }

    // pub fn disperse_encode_rows(&self, data_uncoded: &Vec<Vec<Fr>>) -> Vec<Vec<Fr>> {
    //     let mut data_coded = Vec::new();
    //
    //     let timer_outer = start_timer!(|| "Encoding rows");
    //     for j in 0..self.L {
    //         let timer_inner = start_timer!(|| format!("Row {}", j));
    //
    //         let poly_poly = DensePolynomial::<Fr>::from_coefficients_slice(&data_uncoded[j]);
    //         let poly_evals = poly_poly.evaluate_over_domain(self.domain_encoding);
    //         data_coded.push(poly_evals.evals);
    //
    //         end_timer!(timer_inner);
    //     }
    //     end_timer!(timer_outer);
    //
    //     data_coded
    // }

    // pub fn disperse_encode_rows_systematic(
    //     &self,
    //     data_uncoded: &Vec<Vec<Fr>>,
    // ) -> Vec<Vec<Fr>> {
    //     let mut data_coded = Vec::new();
    //
    //     let timer_outer = start_timer!(|| "Encoding rows");
    //     let domain_uncoded: GeneralEvaluationDomain<Fr> =
    //         ark_poly::domain::EvaluationDomain::<Fr>::new(self.k).unwrap();
    //     for j in 0..self.L {
    //         let timer_inner = start_timer!(|| format!("Row {}", j));
    //         //assert expected length of source data
    //         //assert_eq!(data_uncoded[j].len(), self.k);
    //
    //         let mut poly_evals = Evaluations::from_vec_and_domain(
    //             data_uncoded[j].iter().copied().collect(),
    //             domain_uncoded,
    //         );
    //         let poly_poly = poly_evals.interpolate();
    //         //assert expected degree of interpolated polynomial
    //         //assert_eq!(poly_poly.degree(), self.k-1);
    //
    //         poly_evals = poly_poly.evaluate_over_domain(self.domain_encoding);
    //
    //         data_coded.push(poly_evals.evals);
    //
    //         //assert expected length of erasure coded data
    //         assert_eq!(data_coded[j].len(), self.n);
    //
    //         end_timer!(timer_inner);
    //     }
    //     end_timer!(timer_outer);
    //
    //     //assert uncoded data matches the even indices of coded
    //     // assert_eq!(data_coded[0][0], data_uncoded[0][0]);
    //     //assert_eq!(data_coded[0][2], data_uncoded[0][1]);
    //
    //     data_coded
    // }

    pub fn disperse_encode_rows_lagrange(&self, data_uncoded: &Vec<Vec<Fr>>) -> Vec<Vec<Fr>> {
        let mut data_coded = Vec::with_capacity(data_uncoded.len());

        for row in data_uncoded {
            let mut poly_evals = Vec::with_capacity(self.coded_chunks);
            // Extend with source data first
            poly_evals.extend(row);

            // Then add erasure coded data
            for idx in (0..self.coded_chunks).skip(self.uncoded_chunks) {
                let mut eval = Fr::zero();
                for j in 0..self.uncoded_chunks {
                    let coef = self.evaluate_barycentric(j, idx);
                    eval += row[j] * coef;
                }
                poly_evals.push(eval);
            }

            //assert expected length of erasure coded data
            assert_eq!(poly_evals.len(), self.coded_chunks);

            data_coded.push(poly_evals);
        }

        //assert uncoded data matches the even indices of coded
        // assert_eq!(data_coded[0][0], data_uncoded[0][0]);
        //assert_eq!(data_coded[0][2], data_uncoded[0][1]);

        data_coded
    }

    // pub fn disperse_verify_chunks(
    //     &self,
    //     column_commitments: &Vec<G1Affine>,
    //     data_coded: &Vec<Vec<Fr>>,
    // ) -> bool {
    //     let timer_outer = start_timer!(|| "Checking coded columns");
    //     for i in 0..self.n {
    //         let timer_inner = start_timer!(|| format!("Column {}", i));
    //
    //         let commitment = self.commit_column(&data_coded, i);
    //         let commitment_check = self.encode_commitments(&column_commitments, i);
    //
    //         if commitment != commitment_check {
    //             return false;
    //         }
    //
    //         end_timer!(timer_inner);
    //     }
    //     end_timer!(timer_outer);
    //
    //     true
    // }
    pub fn disperse_verify_chunks_systematic(
        &self,
        source_column_commitments: &Vec<G1Affine>,
        data_coded: &Vec<Vec<Fr>>,
    ) -> bool {
        let timer_outer = start_timer!(|| "Checking coded columns");
        for i in 0..self.coded_chunks {
            let timer_inner = start_timer!(|| format!("Column {}", i));

            let commitment = self.commit_column(&data_coded, i);
            let commitment_interp =
                self.encode_commitments_systematic(&source_column_commitments, i);

            if commitment != commitment_interp {
                return false;
            }

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        true
    }

    pub fn retrieve_download_chunks(
        &self,
        data_coded: &Vec<Vec<Fr>>,
        idxs_download_nodes: &Vec<usize>,
    ) -> Vec<Vec<Fr>> {
        let timer = start_timer!(|| "Downloading chunks");
        let data_coded_downloaded = data_coded
            .iter()
            .map(|r| idxs_download_nodes.iter().map(|&i| r[i]).collect())
            .collect();
        end_timer!(timer);

        data_coded_downloaded
    }

    pub fn retrieve_verify_chunks(
        &self,
        column_commitments: &Vec<G1Affine>,
        data_coded_downloaded: &Vec<Vec<Fr>>,
        idxs_download_nodes: &Vec<usize>,
    ) -> bool {
        let timer_all = start_timer!(|| "Verifying downloaded chunks");

        let timer_encode_commitments =
            start_timer!(|| "'Encoding' of column commitments to coded chunk commitments");
        let mut column_commitments_projective: Vec<G1Projective> = column_commitments
            .iter()
            .map(|h| h.clone().into())
            .collect();
        column_commitments_projective.resize(self.domain_encoding.size(), G1Projective::zero());
        self.domain_encoding
            .fft_in_place(&mut column_commitments_projective);
        let coded_chunk_commitments_affine: Vec<G1Affine> = column_commitments_projective
            .iter()
            .map(|h| h.into_affine())
            .collect();
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

    pub fn retrieve_prepare_decoding(&self, idxs_download_nodes: &Vec<usize>) -> Matrix<Fr> {
        assert!(idxs_download_nodes.len() == self.uncoded_chunks);

        let mut matrix = Vec::new();
        for i in 0..self.uncoded_chunks {
            let i_in_field = Fr::from_le_bytes_mod_order(&i.to_le_bytes());
            matrix.push(
                idxs_download_nodes
                    .iter()
                    .map(|&j| self.domain_encoding.element(j).pow(i_in_field.into_repr()))
                    .collect(),
            );
        }

        Matrix::from_nested_vec(self.uncoded_chunks, self.uncoded_chunks, matrix).invert()
    }

    pub fn retrieve_decode_rows(
        &self,
        data_coded_downloaded: &Vec<Vec<Fr>>,
        decoder_aux: &Matrix<Fr>,
    ) -> Vec<Vec<Fr>> {
        assert!(decoder_aux.height() == decoder_aux.width());
        let mut data_decoded = Vec::new();

        let timer_outer = start_timer!(|| "Decoding rows");
        for j in 0..self.chunk_length {
            let timer_inner = start_timer!(|| format!("Row {}", j));

            assert!(data_coded_downloaded[j].len() == decoder_aux.height());
            data_decoded.push(
                (0..decoder_aux.width())
                    .map(|col| {
                        (0..decoder_aux.height())
                            .map(|row| decoder_aux.get(row, col) * data_coded_downloaded[j][row])
                            .sum::<Fr>()
                    })
                    .collect(),
            );

            end_timer!(timer_inner);
        }
        end_timer!(timer_outer);

        data_decoded
    }

    pub fn sampling_open_entry(
        &self,
        column_commitments: &Vec<G1Affine>,
        data_uncoded: &Vec<Vec<Fr>>,
        row: usize,
        col: usize,
    ) -> (Fr, usize, usize, Vec<G1Affine>, Proof<Bls12_381>) {
        let timer = start_timer!(|| "Poly evaluations and interpolation");
        let poly_evals = Evaluations::from_vec_and_domain(
            data_uncoded.iter().map(|r| r[col]).collect(),
            self.domain_polycommit,
        );
        let poly_poly = poly_evals.interpolate();
        end_timer!(timer);

        let timer = start_timer!(|| "KZG proof");
        // let proof = KZG10::open(&self.kzg10_ck, &poly_poly, self.domain_polycommit.element(row), None).unwrap();
        // Unfortunately, KZG10::open() is pub(crate) only, so inline ... >>>
        let point = self.domain_polycommit.element(row);
        assert!(poly_poly.degree() + 1 <= self.kzg10_ck.size());
        let divisor = DensePolynomial::<Fr>::from_coefficients_vec(vec![-point, Fr::one()]);
        let witness_polynomial = &poly_poly / &divisor;
        assert!(witness_polynomial.degree() + 1 <= self.kzg10_ck.size());
        let proof = Self::unwrap_commitment(
            KZG10::commit(&self.kzg10_ck, &witness_polynomial, None, None).unwrap(),
        );
        let proof = Proof {
            w: proof,
            random_v: None,
        };
        // <<< ... end of inline!
        end_timer!(timer);

        (
            data_uncoded[row][col],
            row,
            col,
            column_commitments.clone(),
            proof,
        )
    }

    pub fn sampling_verify_entry(
        &self,
        (value, row, col, column_commitments, proof): (
            Fr,
            usize,
            usize,
            Vec<G1Affine>,
            Proof<Bls12_381>,
        ),
    ) -> bool {
        let timer = start_timer!(|| "KZG check");
        let commitment = Self::wrap_commitment(column_commitments[col]).0;
        let point = self.domain_polycommit.element(row);
        let ret_val = KZG10::<Bls12_381, DensePolynomial<Fr>>::check(
            &self.kzg10_vk,
            &commitment,
            point,
            value,
            &proof,
        )
        .unwrap();
        end_timer!(timer);

        ret_val
    }
}
