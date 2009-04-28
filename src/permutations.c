/*!
	Generates the next permutation of the vector v of length n.
	@return 1, if there are no more permutations to be generated
	@return 0, otherwise
*/
int next(int v[], int n) {
	/* P2 */
	/* Find the largest i */
	int i = n - 2;
	int jj;
	while ((i >= 0) && (v[i] > v[i + 1])){
		--i;
	}
	/* If i is smaller than 0, then there are no more permutations. */
	if (i < 0){
		return 1;
	}
	/* Find the largest element after vi but not larger than vi */
	int k = n - 1;
	while (v[i] > v[k]){
		--k;
	}
	{       /*   SWAP(v[i], v[k])    */ 
		jj   = v[i];
		v[i] = v[k];
		v[k] = jj;
	}
	/* Swap the last n - i elements. */
	int j;
	k = 0;
	for (j = i + 1; j < (n + i) / 2 + 1; ++j, ++k)
	{/*		SWAP(v[j], v[n - k - 1])    */
		jj = v[j];
		v[j] = v[n - k - 1];
		v[n-k-1] = jj;
	}
	return 0;
}

int flail(int *n, int *v, int *out) {

	/* The initial permutation is 1 2 3 ...*/
	int i,j;
	for (i = 0; i < *n; ++i){
		v[i] = i + 1;
		out[i] = v[i];
	}

	int done = 1;
	i = 0;
	do {
		if (!(done = next(v, *n))){
			i += *n;
			for(j=0; j < *n ; j++){
				out[i+j] = v[j];
			}
		}
	} while (!done);
	return 0;
} 
 
