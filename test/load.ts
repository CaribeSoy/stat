'use strict';

import * as mocha from 'mocha';
describe('Different ways to load util', () => {
	it('Should load the matrix part of the library', () => {
		require('../src/matrix');
	});

	it('Should load the array part of the library', () => {
		require('../src/array');
	});
});