static THRESHOLD: i32 = 30;
static DATA_LENGTH: i32 = 5;

pub fn eval(message: &str) -> Vec<u8> {
    message.as_bytes().to_vec()
}
pub fn erasure_coding(message: &str) -> Vec<(f64, f64)> {
    let message_bytes = eval(message);
    let polynomial_coeffs = message_bytes
        .clone()
        .into_iter()
        .map(|f| f as f64)
        .collect::<Vec<f64>>();
    let mut code_word: Vec<f64> = vec![];
    for i in 0..polynomial_coeffs.len() {
        code_word.push(poly_eval(&polynomial_coeffs, (i + 1) as f64));
    }

    let parity_length = ((THRESHOLD as f64) / 100.0 * code_word.len() as f64).ceil() as usize;

    for i in 0..parity_length {
        code_word.push(poly_eval(
            &polynomial_coeffs,
            (message_bytes.clone().len() + (i + 1)) as f64,
        ))
    }

    let mut result: Vec<(f64, f64)> = vec![];
    for i in 0..code_word.len() {
        result.push(((i + 1) as f64, code_word[i]));
    }
    result
}

pub fn lagrange_interpolation(points: Vec<(f64, f64)>, size: usize) -> Vec<f64> {
    // println!("{:?}", points);
    let n = points.len();
    let mut coeffs = vec![0.0; size]; // Resultant polynomial coefficients

    for i in 0..n {
        let (x_i, y_i) = points[i];
        let mut basis_coeffs = vec![1.0]; // Coefficients for L_i(x)

        for j in 0..n {
            if i != j {
                let (x_j, _) = points[j];

                // Update basis polynomial for L_i(x)
                let mut new_basis = vec![0.0; basis_coeffs.len() + 1];
                for k in 0..basis_coeffs.len() {
                    new_basis[k] -= basis_coeffs[k] * x_j; // -x_j
                    new_basis[k + 1] += basis_coeffs[k]; // +1 * x
                }
                basis_coeffs = new_basis;

                // Scale by the denominator (x_i - x_j)
                let scale = 1.0 / (x_i - x_j);
                for k in 0..basis_coeffs.len() {
                    basis_coeffs[k] *= scale;
                }
            }
        }

        // Add y_i * L_i(x) to the result
        for k in 0..basis_coeffs.len() {
            coeffs[k] += y_i * basis_coeffs[k];
        }
    }

    coeffs.iter().map(|coeff| coeff.abs().round()).collect()
}

pub fn poly_eval(coeffs: &Vec<f64>, point: f64) -> f64 {
    let mut result = 0.0;

    for i in 0..coeffs.len() {
        result += coeffs[i] * point.powi(i.try_into().unwrap());
    }
    result
}

pub fn recover(coeffs: &Vec<(f64, f64)>, point: f64) -> f64 {
    let coeff_y: Vec<_> = coeffs.into_iter().map(|v| v.1).collect();

    let mut result = 0.0;

    for i in 0..coeff_y.len() {
        result += coeff_y[i] * point.powi(i.try_into().unwrap());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = erasure_coding("Hello");

        println!("{:?}", result);
    }

    #[test]
    fn erasure_coding_works_1() {
        let result = erasure_coding("Hello");
        println!("Original evaluations : {:?}", result);
        let malformed_message = vec![
            (1.0, 500.0),
            (3.0, 13254.0),
            (5.0, 86152.0),
            (6.0, 171750.0),
            (7.0, 309626.0),
            (8.0, 517744.0),
        ];
        let recovered_polynomial = lagrange_interpolation(malformed_message, result.len());

        println!("Recovered:{:?}", recovered_polynomial);
        assert_eq!(result[1].1, poly_eval(&recovered_polynomial, 2.0));
    }
    #[test]
    fn erasure_coding_works_2() {
        let message = "Rust";
        let result = erasure_coding(message);
        println!("Original evaluations : {:?}", result);
        let malformed_message = vec![(2.0, 1704.0), (3.0, 4600.0), (4.0, 9814.0), (6.0, 29980.0)];
        let recovered_polynomial = lagrange_interpolation(malformed_message, result.len());

        println!("Recovered coeffs: {:?}", recovered_polynomial);
        assert_eq!(result[1].1, poly_eval(&recovered_polynomial, 2.0));
    }
    #[test]
    fn erasure_coding_works_3() {
        let message = "Propagated";
        let result = erasure_coding(message);
        println!("Original coeffs : {:?}", eval(message));
        println!("Original evaluations : {:?}", result);
        let malformed_message = vec![
            (1.0, 1031.0),
            (2.0, 104608.0),
            (3.0, 2992697.0),
            (5.0, 245743675.0),
            (6.0, 1215364616.0),
            (7.0, 4726557293.0),
            (8.0, 15388807072.0),
            (10.0, 111368394320.0),
            (11.0, 259895541601.0),
            (13.0, 1150627048547.0),
        ];
        let recovered_polynomial = lagrange_interpolation(malformed_message, result.len());

        println!("Recovered coeffs: {:?}", recovered_polynomial);
        assert_eq!(result[1].1, poly_eval(&recovered_polynomial, 2.0));
    }

    #[test]
    fn erasure_coding_works_4() {
        let message = "vitaliketh";
        let result = erasure_coding(message);

        println!("Original evaluations : {:?}", result);
        let malformed_message = vec![
            (1.0, 1077.0),
            (2.0, 109376.0),
            (3.0, 3145357.0),
            (4.0, 37101978.0),
            (5.0, 258411293.0),
            (6.0, 1277163892.0),
            (7.0, 4963322181.0),
            (8.0, 16148603582.0),
            (9.0, 45832084453.0),
            (10.0, 116728689768.0),
            (11.0, 272268587357.0),
            (12.0, 590462098786.0),
            (13.0, 1204389071757.0),
        ];
        let recovered_polynomial = lagrange_interpolation(malformed_message, result.len());

        println!("Recovered coeffs: {:?}", recovered_polynomial);
        assert_eq!(result[1].1, poly_eval(&recovered_polynomial, 2.0));
    }
    #[test]
    fn erasure_coding_works_5() {
        let message = "abcdefghi";
        let result = erasure_coding(message);
       
        let malformed_message = vec![
            (1.0, 909.0),
            (2.0, 53153.0),
            (3.0, 1028389.0),
            (4.0, 9145881.0),
            (5.0, 51147437.0),
            (6.0, 211228489.0),
            (7.0, 705067173.0),
            (8.0, 2010526769.0),
            (9.0, 5078840461.0),
            (10.0, 11654320977.0),
            (11.0, 24734871269.0),
            (12.0, 49209805993.0),
        ];
        let recovered_polynomial = lagrange_interpolation(malformed_message, result.len());

        println!("Recovered coeffs: {:?}", recovered_polynomial);
        assert_eq!(result[1].1, poly_eval(&recovered_polynomial, 2.0));
    }
}
