import * as arrayStat from './array';

/**
 *
 *
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
function compareNumbers(a: number, b: number): number {
	return a - b;
}

/**
 * Computes the maximum of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number}
 */
export function max(matrix: number[][]): number {
	let max = -Infinity;
	for (let i = 0; i < matrix.length; i++) {
		for (let j = 0; j < matrix[i].length; j++) {
			if (matrix[i][j] > max) {
				max = matrix[i][j];
			}
		}
	}
	return max;
}

/**
 * Computes the minimum of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number}
 */
export function min(matrix: number[][]): number {
	let min = Infinity;
	for (let i = 0; i < matrix.length; i++) {
		for (let j = 0; j < matrix[i].length; j++) {
			if (matrix[i][j] < min) {
				min = matrix[i][j];
			}
		}
	}
	return min;
}

/**
 * Computes the min and max of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {{min: number, max: number}}
 */
export function minMax(matrix: number[][]): { min: number, max: number } {
	let min = Infinity;
	let max = -Infinity;
	for (let i = 0; i < matrix.length; i++) {
		for (let j = 0; j < matrix[i].length; j++) {
			if (matrix[i][j] < min) {
				min = matrix[i][j];
			}
			if (matrix[i][j] > max) {
				max = matrix[i][j];
			}
		}
	}
	return {
		min: min,
		max: max
	};
}

/**
 * Computes the entropy of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [eps=0]
 * @returns {number}
 */
export function entropy(matrix: number[][], eps = 0): number {
	let sum = 0,
		l1 = matrix.length,
		l2 = matrix[0].length;
	for (let i = 0; i < l1; i++) {
		for (let j = 0; j < l2; j++) {
			sum += matrix[i][j] * Math.log(matrix[i][j] + eps);
		}
	}
	return -sum;
}

/**
 * Computes the mean of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [dimension=0]
 * @returns {number[]}
 */
export function mean(matrix: number[][], dimension = 0): number[] {
	let rows = matrix.length,
		cols = matrix[0].length,
		theMean: number[], N: number, i: number, j: number;

	if (dimension === -1) {
		theMean = [0];
		N = rows * cols;
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				theMean[0] += matrix[i][j];
			}
		}
		theMean[0] /= N;
	} else if (dimension === 0) {
		theMean = new Array(cols);
		N = rows;
		for (j = 0; j < cols; j++) {
			theMean[j] = 0;
			for (i = 0; i < rows; i++) {
				theMean[j] += matrix[i][j];
			}
			theMean[j] /= N;
		}
	} else if (dimension === 1) {
		theMean = new Array(rows);
		N = cols;
		for (j = 0; j < rows; j++) {
			theMean[j] = 0;
			for (i = 0; i < cols; i++) {
				theMean[j] += matrix[j][i];
			}
			theMean[j] /= N;
		}
	} else {
		throw new Error('Invalid dimension');
	}
	return theMean;
}

/**
 * Computes the sum of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [dimension=0]
 * @returns {number[]}
 */
export function sum(matrix: number[][], dimension = 0): number[] {
	let rows = matrix.length,
		cols = matrix[0].length,
		theSum: number[];

	if (dimension === -1) {
		theSum = [0];
		for (let i = 0; i < rows; i++) {
			for (let j = 0; j < cols; j++) {
				theSum[0] += matrix[i][j];
			}
		}
	} else if (dimension === 0) {
		theSum = new Array(cols);
		for (let j = 0; j < cols; j++) {
			theSum[j] = 0;
			for (let i = 0; i < rows; i++) {
				theSum[j] += matrix[i][j];
			}
		}
	} else if (dimension === 1) {
		theSum = new Array(rows);
		for (let j = 0; j < rows; j++) {
			theSum[j] = 0;
			for (let i = 0; i < cols; i++) {
				theSum[j] += matrix[j][i];
			}
		}
	} else {
		throw new Error('Invalid dimension');
	}
	return theSum;
}

/**
 * Computes the product of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [dimension=0]
 * @returns {number[]}
 */
export function product(matrix: number[][], dimension = 0): number[] {
	let rows = matrix.length,
		cols = matrix[0].length,
		theProduct: number[];

	if (dimension === -1) {
		theProduct = [1];
		for (let i = 0; i < rows; i++) {
			for (let j = 0; j < cols; j++) {
				theProduct[0] *= matrix[i][j];
			}
		}
	} else if (dimension === 0) {
		theProduct = new Array(cols);
		for (let j = 0; j < cols; j++) {
			theProduct[j] = 1;
			for (let i = 0; i < rows; i++) {
				theProduct[j] *= matrix[i][j];
			}
		}
	} else if (dimension === 1) {
		theProduct = new Array(rows);
		for (let j = 0; j < rows; j++) {
			theProduct[j] = 1;
			for (let i = 0; i < cols; i++) {
				theProduct[j] *= matrix[j][i];
			}
		}
	} else {
		throw new Error('Invalid dimension');
	}
	return theProduct;
}

/**
 * Computes the standard deviation of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {boolean} [unbiased]
 * @param {number[]} [means]
 * @returns {number[]}
 */
export function standardDeviation(matrix: number[][], unbiased?: boolean, means?: number[]): number[] {
	let vari = variance(matrix, means, unbiased),
		l = vari.length;
	for (let i = 0; i < l; i++) {
		vari[i] = Math.sqrt(vari[i]);
	}
	return vari;
}

/**
 * Computes the variance of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} [means=mean(matrix)]
 * @param {boolean} [unbiased=true]
 * @returns {number[]}
 */
export function variance(matrix: number[][], means: number[] = mean(matrix), unbiased = true): number[] {
	let rows = matrix.length;
	if (rows === 0) {
		return [];
	}
	let cols = matrix[0].length;
	let vari: number[] = new Array(cols);

	for (let j = 0; j < cols; j++) {
		let sum1 = 0,
			sum2 = 0;
		//x = 0;
		for (let i = 0; i < rows; i++) {
			let x = matrix[i][j] - means[j];
			sum1 += x;
			sum2 += x * x;
		}
		if (unbiased) {
			vari[j] = (sum2 - ((sum1 * sum1) / rows)) / (rows - 1);
		} else {
			vari[j] = (sum2 - ((sum1 * sum1) / rows)) / rows;
		}
	}
	return vari;
}

/**
 * Computes the median of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number[]}
 */
export function median(matrix: number[][]): number[] {
	let rows = matrix.length, cols = matrix[0].length;
	let medians: number[] = new Array(cols);

	for (let i = 0; i < cols; i++) {
		let data = new Array(rows);
		for (let j = 0; j < rows; j++) {
			data[j] = matrix[j][i];
		}
		data.sort(compareNumbers);
		let N = data.length;
		if (N % 2 === 0) {
			medians[i] = (data[N / 2] + data[(N / 2) - 1]) * 0.5;
		} else {
			medians[i] = data[Math.floor(N / 2)];
		}
	}
	return medians;
}

/**
 * Computes the mode of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number[]}
 */
export function mode(matrix: number[][]): number[] {
	let rows = matrix.length,
		cols = matrix[0].length,
		modes: number[] = new Array(cols);

	for (let i = 0; i < cols; i++) {
		let itemCount = new Array(rows);
		for (let k = 0; k < rows; k++) {
			itemCount[k] = 0;
		}
		let itemArray = new Array(rows);
		let count = 0;

		for (let j = 0; j < rows; j++) {
			let index = itemArray.indexOf(matrix[j][i]);
			if (index >= 0) {
				itemCount[index]++;
			} else {
				itemArray[count] = matrix[j][i];
				itemCount[count] = 1;
				count++;
			}
		}

		let maxValue = 0, maxIndex = 0;
		for (let j = 0; j < count; j++) {
			if (itemCount[j] > maxValue) {
				maxValue = itemCount[j];
				maxIndex = j;
			}
		}

		modes[i] = itemArray[maxIndex];
	}
	return modes;
}

/**
 * Computes the skewness of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {boolean} [unbiased=true]
 * @returns {number[]}
 */
export function skewness(matrix: number[][], unbiased = true): number[] {
	let means = mean(matrix);
	let n = matrix.length, l = means.length;
	let skew: number[] = new Array(l);

	for (let j = 0; j < l; j++) {
		let s2 = 0, s3 = 0;
		for (let i = 0; i < n; i++) {
			let dev = matrix[i][j] - means[j];
			s2 += dev * dev;
			s3 += dev * dev * dev;
		}

		let m2 = s2 / n;
		let m3 = s3 / n;
		let g = m3 / Math.pow(m2, 3 / 2);

		if (unbiased) {
			let a = Math.sqrt(n * (n - 1));
			let b = n - 2;
			skew[j] = (a / b) * g;
		} else {
			skew[j] = g;
		}
	}
	return skew;
}

/**
 * Computes the sample kurtosis of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {boolean} [unbiased=true]
 * @returns {number[]}
 */
export function kurtosis(matrix: number[][], unbiased = true): number[] {
	let means = mean(matrix);
	let n = matrix.length, m = matrix[0].length;
	let kurt: number[] = new Array(m);

	for (let j = 0; j < m; j++) {
		let s2 = 0, s4 = 0;
		for (let i = 0; i < n; i++) {
			let dev = matrix[i][j] - means[j];
			s2 += dev * dev;
			s4 += dev * dev * dev * dev;
		}
		let m2 = s2 / n;
		let m4 = s4 / n;

		if (unbiased) {
			let v = s2 / (n - 1);
			let a = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3));
			let b = s4 / (v * v);
			let c = ((n - 1) * (n - 1)) / ((n - 2) * (n - 3));
			kurt[j] = a * b - 3 * c;
		} else {
			kurt[j] = m4 / (m2 * m2) - 3;
		}
	}
	return kurt;
}

/**
 * Computes the standard error of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number[]}
 */
export function standardError(matrix: number[][]): number[] {
	let samples = matrix.length;
	let standardDeviations = standardDeviation(matrix);
	let l = standardDeviations.length;
	let standardErrors: number[] = new Array(l);
	let sqrtN = Math.sqrt(samples);

	for (let i = 0; i < l; i++) {
		standardErrors[i] = standardDeviations[i] / sqrtN;
	}
	return standardErrors;
}

/**
 * Computes the covariance of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [dimension]
 * @returns {number[][]}
 */
export function covariance(matrix: number[][], dimension?: number): number[][] {
	return scatter(matrix, undefined, dimension);
}

/**
 *
 *
 * @export
 * @param {number[][]} matrix
 * @param {number} [divisor=matrix.length - 1]
 * @param {number} [dimension=0]
 * @returns {number[][]}
 */
export function scatter(matrix: number[][], divisor: number = matrix.length - 1, dimension = 0): number[][] {

	if (typeof (divisor) === 'undefined') {
		if (dimension === 0) {
			divisor = matrix.length - 1;
		} else if (dimension === 1) {
			divisor = matrix[0].length - 1;
		}
	}

	let means = mean(matrix, dimension);
	let rows = matrix.length;
	if (rows === 0) {
		return [[]];
	}
	let cols = matrix[0].length,
		cov: number[][], j: number, s: number;

	if (dimension === 0) {
		cov = new Array(cols);
		for (let i = 0; i < cols; i++) {
			cov[i] = new Array(cols);
		}
		for (let i = 0; i < cols; i++) {
			for (j = i; j < cols; j++) {
				s = 0;
				for (let k = 0; k < rows; k++) {
					s += (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
				}
				s /= divisor;
				cov[i][j] = s;
				cov[j][i] = s;
			}
		}
	} else if (dimension === 1) {
		cov = new Array(rows);
		for (let i = 0; i < rows; i++) {
			cov[i] = new Array(rows);
		}
		for (let i = 0; i < rows; i++) {
			for (j = i; j < rows; j++) {
				s = 0;
				for (let k = 0; k < cols; k++) {
					s += (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
				}
				s /= divisor;
				cov[i][j] = s;
				cov[j][i] = s;
			}
		}
	} else {
		throw new Error('Invalid dimension');
	}

	return cov;
}

/**
 * Computes the correlation of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @returns {number[][]}
 */
export function correlation(matrix: number[][]): number[][] {
	let means = mean(matrix),
		standardDeviations = standardDeviation(matrix, true, means),
		scores = zScores(matrix, means, standardDeviations),
		rows = matrix.length,
		cols = matrix[0].length;

	let cor: number[][] = new Array(cols);
	for (let i = 0; i < cols; i++) {
		cor[i] = new Array(cols);
	}
	for (let i = 0; i < cols; i++) {
		for (let j = i; j < cols; j++) {
			let c = 0;
			for (let k = 0, l = scores.length; k < l; k++) {
				c += scores[k][j] * scores[k][i];
			}
			c /= rows - 1;
			cor[i][j] = c;
			cor[j][i] = c;
		}
	}
	return cor;
}

/**
 * Computes the zScores of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} [means=mean(matrix)]
 * @param {boolean} [standardDeviations=standardDeviation(matrix, true, means)]
 * @returns {number[][]}
 */
export function zScores(matrix: number[][], means: number[] = mean(matrix), standardDeviations = standardDeviation(matrix, true, means)): number[][] {
	return standardize(center(matrix, means, false), standardDeviations, true);
}

/**
 * Computes the center of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} [means=mean(matrix)]
 * @param {boolean} inPlace
 * @returns {number[][]}
 */
export function center(matrix: number[][], means: number[] = mean(matrix), inPlace: boolean): number[][] {
	let result = matrix,
		l = matrix.length;

	if (!inPlace) {
		result = new Array(l);
		for (let i = 0; i < l; i++) {
			result[i] = new Array(matrix[i].length);
		}
	}

	for (let i = 0; i < l; i++) {
		let row = result[i];
		for (let j = 0, jj = row.length; j < jj; j++) {
			row[j] = matrix[i][j] - means[j];
		}
	}
	return result;
}

/**
 * Standardize the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} [standardDeviations=standardDeviation(matrix)]
 * @param {boolean} [inPlace=true]
 * @returns {number[][]}
 */
export function standardize(matrix: number[][], standardDeviations: number[] = standardDeviation(matrix), inPlace = true): number[][] {
	let result = matrix,
		l = matrix.length;

	if (!inPlace) {
		result = new Array(l);
		for (let i = 0; i < l; i++) {
			result[i] = new Array(matrix[i].length);
		}
	}

	for (let i = 0; i < l; i++) {
		let resultRow = result[i];
		let sourceRow = matrix[i];
		for (let j = 0, jj = resultRow.length; j < jj; j++) {
			if (standardDeviations[j] !== 0 && !isNaN(standardDeviations[j])) {
				resultRow[j] = sourceRow[j] / standardDeviations[j];
			}
		}
	}
	return result;
}

/**
 * Computes the weighted variance of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @returns {number[]}
 */
export function weightedVariance(matrix: number[][], weights: number[]): number[] {
	let means = mean(matrix);
	let rows = matrix.length;
	if (rows === 0) {
		return [];
	}
	let cols = matrix[0].length;
	let vari: number[] = new Array(cols);

	for (let j = 0; j < cols; j++) {
		let sum = 0;
		let a = 0, b = 0;

		for (let i = 0; i < rows; i++) {
			let z = matrix[i][j] - means[j];
			let w = weights[i];

			sum += w * (z * z);
			b += w;
			a += w * w;
		}

		vari[j] = sum * (b / (b * b - a));
	}

	return vari;
}

/**
 * Computes the weighted mean of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param {number} [dimension=0]
 * @returns {number[]}
 */
export function weightedMean(matrix: number[][], weights: number[], dimension = 0): number[] {
	let rows = matrix.length;
	if (rows === 0) {
		return [];
	}
	let cols = matrix[0].length,
		means: number[], w: number, row: number[];

	if (dimension === 0) {
		means = new Array(cols);
		for (let i = 0; i < cols; i++) {
			means[i] = 0;
		}
		for (let i = 0; i < rows; i++) {
			row = matrix[i];
			w = weights[i];
			for (let j = 0; j < cols; j++) {
				means[j] += row[j] * w;
			}
		}
	} else if (dimension === 1) {
		means = new Array(rows);
		for (let i = 0; i < rows; i++) {
			means[i] = 0;
		}
		for (let j = 0; j < rows; j++) {
			row = matrix[j];
			w = weights[j];
			for (let i = 0; i < cols; i++) {
				means[j] += row[i] * w;
			}
		}
	} else {
		throw new Error('Invalid dimension');
	}

	let weightSum = arrayStat.sum(weights);
	if (weightSum !== 0) {
		for (let i = 0, ii = means.length; i < ii; i++) {
			means[i] /= weightSum;
		}
	}
	return means;
}

/**
 * Computes the weighted covariance of the given matrix
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param {number} [dimension=0]
 * @param {number[]} [means=weightedMean(matrix, weights, dimension)]
 * @returns {number[][]}
 */
export function weightedCovariance(matrix: number[][], weights: number[], dimension = 0, means: number[] = weightedMean(matrix, weights, dimension)): number[][] {
	let s1 = 0, s2 = 0;
	for (let i = 0, ii = weights.length; i < ii; i++) {
		s1 += weights[i];
		s2 += weights[i] * weights[i];
	}
	let factor = s1 / (s1 * s1 - s2);
	return weightedScatter(matrix, weights, means, factor, dimension);
}

/**
 *
 *
 * @export
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param {number[]} means
 * @param {number} [factor=1]
 * @param {number} [dimension=0]
 * @returns {number[][]}
 */
export function weightedScatter(matrix: number[][], weights: number[], means: number[], factor = 1, dimension = 0): number[][] {
	means = means || weightedMean(matrix, weights, dimension);
	let rows = matrix.length;
	if (rows === 0) {
		return [[]];
	}
	let cols = matrix[0].length,
		cov: number[][];

	if (dimension === 0) {
		cov = new Array(cols);
		for (let i = 0; i < cols; i++) {
			cov[i] = new Array(cols);
		}
		for (let i = 0; i < cols; i++) {
			for (let j = i; j < cols; j++) {
				let s = 0;
				for (let k = 0; k < rows; k++) {
					s += weights[k] * (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
				}
				cov[i][j] = s * factor;
				cov[j][i] = s * factor;
			}
		}
	} else if (dimension === 1) {
		cov = new Array(rows);
		for (let i = 0; i < rows; i++) {
			cov[i] = new Array(rows);
		}
		for (let i = 0; i < rows; i++) {
			for (let j = i; j < rows; j++) {
				let s = 0;
				for (let k = 0; k < cols; k++) {
					s += weights[k] * (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
				}
				cov[i][j] = s * factor;
				cov[j][i] = s * factor;
			}
		}
	} else {
		throw new Error('Invalid dimension');
	}
	return cov;
}