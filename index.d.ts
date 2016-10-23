// Type definitions for ml-stat v1.3.3
// Project: https://github.com/mljs/stat
// Definitions by: Diana Caro <https://github.com/mlts/stat>


declare module 'ml-stat' {
    export import array = require('ml-stat/array');
    export import matrix = require('ml-stat/matrix');
}

declare module 'ml-stat/array' {
    /**
        * Computes the sum of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function sum(values: number[]): number;
    /**
        * Computes the maximum of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function max(values: number[]): number;
    /**
        * Computes the minimum of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function min(values: number[]): number;
    /**
        * Computes the min and max of the given values
        * @param {Array} values
        * @returns {{min: number, max: number}}
        */
    export function minMax(values: number[]): {
            min: number;
            max: number;
    };
    /**
        * Computes the arithmetic mean of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function arithmeticMean(values: number[]): number;
    /**
        * {@link arithmeticMean}
        */
    export var mean: typeof arithmeticMean;
    /**
        * Computes the geometric mean of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function geometricMean(values: number[]): number;
    /**
        * Computes the mean of the log of the given values
        * If the return value is exponentiated, it gives the same result as the
        * geometric mean.
        * @param {Array} values
        * @returns {number}
        */
    export function logMean(values: number[]): number;
    /**
        * Computes the weighted grand mean for a list of means and sample sizes
        * @param {Array} means - Mean values for each set of samples
        * @param {Array} samples - Number of original values for each set of samples
        * @returns {number}
        */
    export function grandMean(means: number[], samples: number[]): number;
    /**
        * Computes the truncated mean of the given values using a given percentage
        * @param {Array} values
        * @param {number} percent - The percentage of values to keep (range: [0,1])
        * @param {boolean} [alreadySorted=false]
        * @returns {number}
        */
    export function truncatedMean(values: number[], percent: number, alreadySorted?: boolean): number;
    /**
        * Computes the harmonic mean of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function harmonicMean(values: number[]): number;
    /**
        * Computes the contra-harmonic mean of the given values
        * @param {Array} values
        * @returns {number}
        */
    export function contraHarmonicMean(values: number[]): number;
    /**
        * Computes the median of the given values
        * @param {Array} values
        * @param {boolean} [alreadySorted=false]
        * @returns {number}
        */
    export function median(values: number[], alreadySorted?: boolean): number;
    /**
        * Computes the variance of the given values
        * @param {Array} values
        * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
        * @returns {number}
        */
    export function variance(values: number[], unbiased?: boolean): number;
    /**
        * Computes the standard deviation of the given values
        * @param {Array} values
        * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
        * @returns {number}
        */
    export function standardDeviation(values: number[], unbiased?: boolean): number;
    export function standardError(values: number[]): number;
    /**
        * IEEE Transactions on biomedical engineering, vol. 52, no. 1, January 2005, p. 76-
        * Calculate the standard deviation via the Median of the absolute deviation
        *  The formula for the standard deviation only holds for Gaussian random variables.
        * @returns {{mean: number, stdev: number}}
        */
    export function robustMeanAndStdev(y: number[]): {
            mean: number;
            stdev: number;
    };
    export function quartiles(values: number[], alreadySorted?: boolean): {
            q1: number;
            q2: number;
            q3: number;
    };
    export function pooledStandardDeviation(samples: number[][], unbiased?: boolean): number;
    export function pooledVariance(samples: number[][], unbiased?: boolean): number;
    export function mode(values: number[]): any;
    export function covariance(vector1: number[], vector2: number[], unbiased?: boolean): number;
    export function skewness(values: number[], unbiased?: boolean): number;
    export function kurtosis(values: number[], unbiased?: boolean): number;
    export function entropy(values: number[], eps: number): number;
    export function weightedMean(values: number[], weights: number[]): number;
    export function weightedStandardDeviation(values: number[], weights: number[]): number;
    export function weightedVariance(values: number[], weights: number[]): number;
    export function center(values: number[], inPlace?: boolean): void;
    export function standardize(values: number[], standardDev?: number, inPlace?: boolean): any[];
    export function cumulativeSum(array: number[]): any[];
}

declare module 'ml-stat/matrix' {
    export function max(matrix: number[][]): number;
    export function min(matrix: number[][]): number;
    export function minMax(matrix: number[][]): {
        min: number;
        max: number;
    };
    export function entropy(matrix: number[][], eps?: number): number;
    export function mean(matrix: number[][], dimension?: number): number[];
    export function sum(matrix: number[][], dimension?: number): number[];
    export function product(matrix: number[][], dimension?: number): number[];
    export function standardDeviation(matrix: number[][], means?: number[], unbiased?: boolean): number[];
    export function variance(matrix: number[][], means: number[], unbiased?: boolean): number[];
    export function median(matrix: number[][]): number[];
    export function mode(matrix: number[][]): any[];
    export function skewness(matrix: number[][], unbiased?: boolean): number[];
    export function kurtosis(matrix: number[][], unbiased?: boolean): any[];
    export function standardError(matrix: number[][]): any[];
    export function covariance(matrix: number[][], dimension: number): number[][];
    export function scatter(matrix: number[][], divisor: number, dimension?: number): number[][];
    export function correlation(matrix: number[][]): any[];
    export function zScores(matrix: number[][], means: number[], standardDeviations?: number[]): number[][];
    export function center(matrix: number[][], means: number[], inPlace: boolean): number[][];
    export function standardize(matrix: number[][], standardDeviations: number[], inPlace: boolean): number[][];
    export function weightedVariance(matrix: number[][], weights: number[]): number[];
    export function weightedMean(matrix: number[][], weights: number[], dimension?: number): number[];
    export function weightedCovariance(matrix: number[][], weights: number[], means: number[], dimension: number): number[][];
    export function weightedScatter(matrix: number[][], weights: number[], means: number[], factor: number, dimension: number): number[][];
}

