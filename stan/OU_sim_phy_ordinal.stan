functions { 
matrix kronecker_prod(matrix A, matrix B) { // returns the Kronecker Product
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m)
    for (j in 1:n)
      for (k in 1:p)
        for (l in 1:q)
      C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];
  return C;
}

matrix A_dt(matrix A, real t) {  // expected auto and cross effects over a discrete time t
 return( matrix_exp(A * t) );
}

matrix A_sharp(matrix A) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_temp;
  matrix[rows(A),cols(A)] I; // identity matrix
  
  I = diag_matrix(rep_vector(1,rows(A)));

  A_temp = kronecker_prod(A,I) + kronecker_prod(I,A);
  return(A_temp);
}

matrix cov_drift(matrix A, matrix Q, real ts) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_sharp_temp;
  matrix[rows(A) * rows(A), cols(A) * cols(A)] I; // identity matrix
  vector[rows(Q)*cols(Q)] row_Q;
  vector[rows(A)*cols(A)] irow_vec;
  matrix[rows(A),cols(A)] irow_mat;
  
  I = diag_matrix(rep_vector(1,rows(A_sharp_temp)));
  
  A_sharp_temp = A_sharp(A);
  
  // row operation takes elements of a matrix rowwise and puts them into a column vector
  for (i in 1:rows(Q))
  for (j in 1:cols(Q)) {
    row_Q[i + (j-1)*rows(Q)] = Q[j,i];
  }
  
  irow_vec = inverse(A_sharp_temp) * (matrix_exp(A_sharp_temp * ts) - I) * row_Q;
  
  // irow takes elements of a column vector and puts them in a matrix rowwise
  {
  int row_size = rows(A);
  int row_ticker = 1;
  int col_ticker = 0;
  
  for (i in 1:num_elements(irow_vec)) {
    col_ticker += 1;
    if (col_ticker > row_size) {
      row_ticker += 1;
      col_ticker = 1;
    }
    irow_mat[row_ticker,col_ticker] = irow_vec[i];
  }
  }
  return(irow_mat);
}

} // end of functions block


data{
  int N;
  int J;
  int N_seg;
  int node_seq[N_seg];
  int parent[N_seg];
  vector[N_seg] ts;
}

parameters{
vector[2] b; // SDE intercepts
vector[2] eta_anc; // ancestral states
matrix[N_seg - 1,2] z_drift; // stochastic drift, unscaled and uncorrelated

ordered[3] c1; // cutpoints for observed variable 1
ordered[3] c2; // cutpoints for observed variable 2
}

transformed parameters{
matrix[N_seg,2] eta;
matrix[2,2] Q;
matrix[2,2] I; 
matrix[2,2] A;

// selection matrix, col 1 is var 1, col 2 is var 2
A[1,1] = -1;   // autoregressive effect of var 1 on itself
A[2,2] = -0.5; // autoregressive effect of var 2 on itself

A[2,1] = 4;
A[1,2] = 0;

// drift matrix
Q[1,1] = 6;
Q[2,2] = 3;
Q[1,2] = 0;
Q[2,1] = 0;

// identity matrix
I = diag_matrix(rep_vector(1,2));

// setting ancestral states
eta[node_seq[1],1] = eta_anc[1]; // var 1 ancestral state
eta[node_seq[1],2] = eta_anc[2]; // var 2 ancestral state

  for (i in 2:N_seg) {
  matrix[2,2] A_delta; // amount of deterministic change ("selection")
  matrix[2,2] VCV; // variance-covariance matrix of stochastic change ("drift")
  vector[2] drift_seg; // accumulated drift over the segment
  
  A_delta = A_dt(A, ts[i]);
  VCV = cov_drift(A, Q, ts[i]);
  
  // No drift on the interaction, bc its simply a product of var 1 and var 2
  drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );
  
  eta[node_seq[i],] = to_row_vector(
    A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg
  );
  }

} // end transformed parameters block

model{
  b ~ std_normal();
  eta_anc ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  c1 ~ normal(0, 2);
  c2 ~ normal(0, 2);
}

generated quantities {
  int y[N,2];
  
  for (i in 1:N) {
    y[i,1] = ordered_logistic_rng( eta[i,1], c1 );
    y[i,2] = ordered_logistic_rng( eta[i,2], c2 );
  }
}
