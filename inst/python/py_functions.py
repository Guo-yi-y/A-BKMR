import numpy as np
import mpmath
import pandas as pd

from scipy.spatial.distance import cdist



def euclideanDist(x, y):
    return cdist(x, y, 'euclidean')

def sam_py(R, nd, id_ini, num_nn=10, P=-20, Q=20, max_loop=20):
    R = np.array(R)  
    nd = int(nd)
    id = np.arange(R.shape[0])
    

    Dset = np.array(id_ini)
    #Dset = np.random.choice(id, nd, replace=False)
    Cset = np.setdiff1d(id, Dset)
    

    dist_mat = euclideanDist(R[Cset, :], R[Dset, :])
    rs = np.sum(dist_mat**P, axis=1)
    crit_i = np.sum(rs ** (Q / P)) ** (1 / Q)
    

    CRIT = np.full(len(Cset), np.inf)
    
    loop_counter = 0
    crit_old = np.inf
    
    while crit_i != crit_old and loop_counter < max_loop:
        for i in range(nd):
            Dset_i = R[Dset[i], np.newaxis]
            
            partial_newrow = np.sum(euclideanDist(Dset_i, R[np.setdiff1d(Dset, Dset[i]), :]) ** P)
            rs_without_i = rs - dist_mat[:, i] ** P
            

            vec = np.argsort(dist_mat[:, i])[:num_nn]
            
            for j in vec:
                Cset_j = R[Cset[j], np.newaxis]
                newcol = (euclideanDist(Cset_j, R[np.setdiff1d(id, np.append(Dset, Cset[j])), :])) ** P
                
                CRIT[j] = (np.sum((np.delete(rs_without_i, j) + newcol) ** (Q / P)) + 
                           np.sum((euclideanDist(Cset_j, Dset_i) ** P + partial_newrow) ** (Q / P)) ) ** (1 / Q)
            

            best = np.min(CRIT)
            best_spot = Cset[np.argmin(CRIT)]
            crit_old = crit_i
            if best < crit_i:
                crit_i = best
                Dset[i] = best_spot
                Cset = np.setdiff1d(id, Dset)
                dist_mat = euclideanDist(R[Cset, :], R[Dset, :])
                rs = np.sum(dist_mat**P, axis=1)
        
        loop_counter += 1
        #print(Dset+1)
    
    return Dset+1

def sam_py_w(R, nd, id_ini, num_nn=10, P=-20, Q=20, max_loop=20):
    R_ic = pd.DataFrame(R)
    R_ic['id'] = pd.factorize(R_ic.astype(str).agg('-'.join, axis=1))[0]
    R_uni = R_ic.groupby(list(R_ic.columns)).size().reset_index(name='count')
    R_uni = R_uni.sort_values('id')
    R = np.array(R_uni)  
    
    
    id_loc = R_uni.columns.get_loc('id')
    weight_loc = R_uni.columns.get_loc('count')
    re_loc = [id_loc, weight_loc]
    nd = int(nd)
    id = np.array(R[:,id_loc]).astype(int)
    
    p_sam = R[:, weight_loc]/np.sum(R[:, weight_loc])

    Dset = np.array(id_ini)
    #random.seed(111)
    #Dset = np.random.choice(id, nd, replace=False, p = p_sam)
    #Dset = np.arange(nd)
    Cset = np.setdiff1d(id, Dset)
    

    ne_loc = np.array([i for i in range(R.shape[1]) if i not in [weight_loc, id_loc]])
    dist_mat = euclideanDist(R[Cset, :][:, ne_loc], R[Dset, :][:, ne_loc])
    weight_mat = np.outer(R[Cset, weight_loc], R[Dset, weight_loc])
    rs = np.sum(dist_mat**P*weight_mat, axis=1)
    crit_i = np.sum(rs ** (Q / P)) ** (1 / Q)
    

    CRIT = np.full(len(Cset), np.inf)
    
    loop_counter = 0
    crit_old = np.inf
    
    while crit_i != crit_old and loop_counter < max_loop:
        for i in range(nd):
            Dset_i = R[Dset[i], np.newaxis][:, ne_loc]
            weight_i = R[Dset[i], weight_loc]
            
            partial_newrow = np.sum(euclideanDist(Dset_i, R[np.setdiff1d(Dset, Dset[i]), :][:, ne_loc]) ** P *np.outer(weight_i, R[np.setdiff1d(Dset, Dset[i]), weight_loc]))
            rs_without_i = rs - dist_mat[:, i] ** P*np.outer(weight_i, R[Cset, weight_loc])
            

            vec = np.argsort(dist_mat[:, i]*weight_mat[:, i])[:num_nn]
            
            for j in vec:
                Cset_j = R[Cset[j], np.newaxis][:, ne_loc]
                weight_j = R[Cset[j], weight_loc]
                newcol = (euclideanDist(Cset_j, R[np.setdiff1d(id, np.append(Dset, Cset[j])), :][:, ne_loc])) ** P*np.outer(weight_j, R[np.setdiff1d(id, np.append(Dset, Cset[j])), weight_loc]  )
                
                CRIT[j] = (np.sum((np.delete(rs_without_i, j) + newcol) ** (Q / P)) + 
                           np.sum((euclideanDist(Cset_j, Dset_i) ** P*np.outer(weight_i, weight_j) + partial_newrow) ** (Q / P)) ) ** (1 / Q)
            

            best = np.min(CRIT)
            best_spot = Cset[np.argmin(CRIT)]
            crit_old = crit_i
            if best < crit_i:
                crit_i = best
                Dset[i] = best_spot
                Cset = np.setdiff1d(id, Dset)
                dist_mat = euclideanDist(R[Cset, :][:, ne_loc], R[Dset, :][:, ne_loc])
                weight_mat = np.outer(R[Cset, weight_loc], R[Dset, weight_loc])
                rs = np.sum(dist_mat**P*weight_mat, axis=1)
        
        loop_counter += 1
        #print(Dset+1)
    
    return Dset+1

def exp_neg_mat(matrix):
    """Applies a negative logarithm transformation to the matrix."""
    return np.exp(-matrix)


def tcrossprod_py(matrix):
    """
    Compute the cross-product of the transpose of a matrix with the matrix itself.
    
    Equivalent to `tcrossprod(matrix)` in R.
    """
    return np.dot(matrix, matrix.T)

def comp_XVinv(X, lambda_1, K10, Rinv):
    """
    Performs the matrix operation t(X) - lambda[1] * t(X) @ t(K10) @ Rinv @ K10.
    """
    X_t = X.T
    K10_t = K10.T
    intermediate = X_t @ K10_t @ Rinv @ K10
    result = X_t - lambda_1 * intermediate
    return result



def compute_muVinv(mu, lambda_val, K10, Rinv):

    
    muVinv = np.dot(mu.T, mu) - lambda_val * np.dot(np.dot(np.dot(mu.T, K10.T), Rinv), np.dot(K10, mu))
    return muVinv
  
  
def compute_Vinv(lambda_val, K10, Rinv):

    
    Vinv = lambda_val * np.dot(np.dot(K10.T, Rinv), K10)
    return Vinv
  
def compute_stainv(mu, lambda_val, K10, Rinv):
  stain = mu.T - lambda_val*np.dot(np.dot(np.dot(mu.T, K10.T), Rinv), K10)
  
  return stain
