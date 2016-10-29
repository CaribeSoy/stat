// Type definitions for ml-stat v1.3.3
// Project: https://github.com/mljs/stat
// Definitions by: Diana Caro <https://github.com/CaribeSoy>

declare module "ml-stat" {

	namespace array {
		/**
		 * Computes the sum of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function sum(values: number[]): number;

		/**
		 * Computes the maximum of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function max(values: number[]): number;

		/**
		 * Computes the minimum of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function min(values: number[]): number;

		/**
		 * Computes the min and max of the given values
		 * @param {Array} values
		 * @returns {{min: number, max: number}}
		 */
		function minMax(values: number[]): {
			min: number;
			max: number;
		};

		/**
		 * Computes the arithmetic mean of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function arithmeticMean(values: number[]): number;

		/**
		 * {@link arithmeticMean}
		 */
		var mean: typeof arithmeticMean;

		/**
		 * Computes the geometric mean of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function geometricMean(values: number[]): number;

		/**
		 * Computes the mean of the log of the given values
		 * If the return value is exponentiated, it gives the same result as the
		 * geometric mean.
		 * @param {Array} values
		 * @returns {number}
		 */
		function logMean(values: number[]): number;

		/**
		 * Computes the weighted grand mean for a list of means and sample sizes
		 * @param {Array} means - Mean values for each set of samples
		 * @param {Array} samples - Number of original values for each set of samples
		 * @returns {number}
		 */
		function grandMean(means: number[], samples: number[]): number;

		/**
		 * Computes the truncated mean of the given values using a given percentage
		 * @param {Array} values
		 * @param {number} percent - The percentage of values to keep (range: [0,1])
		 * @param {boolean} [alreadySorted=false]
		 * @returns {number}
		 */
		function truncatedMean(values: number[], percent: number, alreadySorted?: boolean): number;

		/**
		 * Computes the harmonic mean of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function harmonicMean(values: number[]): number;

		/**
		 * Computes the contra-harmonic mean of the given values
		 * @param {Array} values
		 * @returns {number}
		 */
		function contraHarmonicMean(values: number[]): number;

		/**
		 * Computes the median of the given values
		 * @param {Array} values
		 * @param {boolean} [alreadySorted=false]
		 * @returns {number}
		 */
		function median(values: number[], alreadySorted?: boolean): number;
    
		/**
		 * Computes the variance of the given values
		 * @param {Array} values
		 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
		 * @returns {number}
		 */
		function variance(values: number[], unbiased?: boolean): number;

		/**
		 * Computes the standard deviation of the given values
		 * @param {Array} values
		 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
		 * @returns {number}
		 */
		function standardDeviation(values: number[], unbiased?: boolean): number;

		function standardError(values: number[]): number;

		/**
		 * IEEE Transactions on biomedical engineering, vol. 52, no. 1, January 2005, p. 76-
		 * Calculate the standard deviation via the Median of the absolute deviation
		 *  The formula for the standard deviation only holds for Gaussian random variables.
		 * @returns {{mean: number, stdev: number}}
		 */
		function robustMeanAndStdev(y: number[]): {
			mean: number;
			stdev: number;
		};

		function quartiles(values: number[], alreadySorted?: boolean): {
			q1: number;
			q2: number;
			q3: number;
		};

		function pooledStandardDeviation(samples: number[][], unbiased?: boolean): number;

		function pooledVariance(samples: number[][], unbiased?: boolean): number;

		/**
		 *
		 * @param {number[]} values
		 * @returns
		 */
		function mode(values: number[]): any;

		/**
		 *
		 * @param {number[]} vector1
		 * @param {number[]} vector2
		 * @param unbiased = true
		 * @returns
		 */
		function covariance(vector1: number[], vector2: number[], unbiased?: boolean): number;

		/**
		 *
		 * @param {number[]} values
		 * @param unbiased = true
		 * @returns
		 */
		function skewness(values: number[], unbiased?: boolean): number;

		/**
		 *
		 * @param {number[]} values
		 * @param unbiased = true
		 * @returns
		 */
		function kurtosis(values: number[], unbiased?: boolean): number;

		/**
		 *
		 * @param {number[]} values
		 * @param eps = 0
		 * @returns
		 */
		function entropy(values: number[], eps?: number): number;

		/**
		 *
		 * @param {number[]} values
		 * @param {number[]} weights
		 * @returns
		 */
		function weightedMean(values: number[], weights: number[]): number;

		/**
		 *
		 * @param {number[]} values
		 * @param {number[]} weights
		 * @returns
		 */
		function weightedStandardDeviation(values: number[], weights: number[]): number;

		/**
		 *
		 * @param {number[]} values
		 * @param {number[]} weights
		 * @returns
		 */
		function weightedVariance(values: number[], weights: number[]): number;

		/**
		 *
		 * @param {number[]} values
		 * @param inPlace = false
		 * @returns
		 */
		function center(values: number[], inPlace?: boolean): number[];

		/**
		 *
		 * @param {number[]} values
		 * @param standardDev = standardDeviation(values)
		 * @param inPlace = false
		 * @returns
		 */
		function standardize(values: number[], standardDev?: number, inPlace?: boolean): number[];

		/**
		 *
		 * @param {number[]} array
		 * @returns
		 */
		function cumulativeSum(array: number[]): number[];
	}

	namespace matrix {
		/**
		 *
		 * @param {number[][]} matrix
		 * @returns
		 */
		function max(matrix: number[][]): number;

		/**
		 *
		 * @param {number[][]} matrix
		 * @returns
		 */
		function min(matrix: number[][]): number;

		/**
		 *
		 * @param {number[][]} matrix
		 */
		function minMax(matrix: number[][]): {
			min: number;
			max: number;
		};

		/**
		 *
		 * @param {number[][]} matrix
		 * @param eps = 0
		 * @returns
		 */
		function entropy(matrix: number[][], eps?: number): number;

		/**
		 *
		 * @param {number[][]} matrix
		 * @param dimension = 0
		 * @returns
		 */
		function mean(matrix: number[][], dimension?: number): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param dimension = 0
		 * @returns
		 */
		function sum(matrix: number[][], dimension?: number): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param dimension = 0
		 * @returns
		 */
		function product(matrix: number[][], dimension?: number): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} means?
		 * @param {boolean} unbiased?
		 * @returns
		 */
		function standardDeviation(matrix: number[][], unbiased?: boolean, means?: number[]): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} means
		 * @param unbiased = true
		 * @returns
		 */
		function variance(matrix: number[][], means?: number[], unbiased?: boolean): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @returns
		 */
		function median(matrix: number[][]): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @returns
		 */
		function mode(matrix: number[][]): number[];

		/**
		 *
		 * @param matrix
		 * @param unbiased
		 */
		function skewness(matrix: number[][], unbiased?: boolean): number[];

		/**
		 *
		 * @param matrix
		 * @param unbiased
		 */
		function kurtosis(matrix: number[][], unbiased?: boolean): number[];

		/**
		 *
		 * @param matrix
		 */
		function standardError(matrix: number[][]): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number} dimension?
		 * @returns
		 */
		function covariance(matrix: number[][], dimension?: number): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number} divisor?
		 * @param dimension = 0
		 * @returns
		 */
		function scatter(matrix: number[][], divisor?: number, dimension?: number): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @returns
		 */
		function correlation(matrix: number[][]): number[][];

		/**
		 *
		 * @param matrix
		 * @param means
		 * @param standardDeviations
		 */
		function zScores(matrix: number[][], means?: number[], standardDeviations?: number[]): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} means
		 * @param {boolean} inPlace
		 * @returns
		 */
		function center(matrix: number[][], means: number[], inPlace: boolean): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param standardDeviations = standardDeviation(matrix)
		 * @param {boolean} inPlace
		 * @returns
		 */
		function standardize(matrix: number[][], standardDeviations: number[], inPlace: boolean): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} weights
		 * @returns
		 */
		function weightedVariance(matrix: number[][], weights: number[]): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} weights
		 * @param dimension = 0
		 * @returns
		 */
		function weightedMean(matrix: number[][], weights: number[], dimension?: number): number[];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} weights
		 * @param dimension = 0
		 * @param means = weightedMean(matrix
		 * @param weights
		 * @param dimension)
		 * @returns
		 */
		function weightedCovariance(matrix: number[][], weights: number[], dimension?: number, means?: number[]): number[][];

		/**
		 *
		 * @param {number[][]} matrix
		 * @param {number[]} weights
		 * @param {number[]} means
		 * @param factor = 1
		 * @param {number} dimension?
		 * @returns
		 */
		function weightedScatter(matrix: number[][], weights: number[], means: number[], factor?: number, dimension?: number): number[][];

	}
}
