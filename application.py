import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
# os.chdir(r"C:\Users\USER\Dropbox\Python exercise\Python code for chapter 3")




def sBIC(y, sep_b, sep_c):    
    
    Big_num = np.inf
    s2_min_UB, s2_min_UBU, s2_min_UBC, s2_min_UBCU = Big_num, Big_num, Big_num, Big_num
    tb_hat_UB, tb_hat_UBU, tb_hat_UBC, tb_hat_UBCU = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

    n=len(y)
    
    fd_y=y-pd.Series(y).shift(1)
    sign_y=(fd_y/np.abs(fd_y)).to_numpy().reshape(-1, 1)

    x_mat=np.ones((n,1))
    # ssrU=(sign_y[1:n,].T @ sign_y[1:n,])/n
    
    

    ssrjT_U, ssrjT_C, ssr1T = np.zeros((n,1)), np.zeros((n,1)), np.zeros((n,1))


    for t1 in range(sep_b, n-1):
        
        ssrjT_U[t1, 0] = (sign_y[t1+1:n,:].T @ sign_y[t1+1:n,:]).item()

        if t1 < n - 2:
            dyls=sign_y[t1+1:n]
            xs=x_mat[t1+1:n,:]
            ixsxs = np.linalg.inv(xs.T @ xs)
            z = xs.T @ dyls
            ssrjT_C[t1, 0]=(dyls.T @ dyls - z.T @ ixsxs @ z).item()


    for t1 in range(sep_b-1, n-1):
        ssr01 = (sign_y[1:t1+1,:].T @ sign_y[1:t1+1,:]).item()

        if y[-1] > y[t1] and t1 < n - 2:
            dyls=sign_y[t1+1:n,:]
            xs=x_mat[t1+1:n,:]
            ixsxs = np.linalg.inv(xs.T @ xs)
            z = xs.T @ dyls
            ssr1T=(dyls.T @ dyls - z.T @ ixsxs @ z).item()         
            s2 = (ssr01 + ssr1T) / n
            if s2 < s2_min_UB: 
                s2_min_UB, tb_hat_UB[0], tb_hat_UB[1] = s2, t1+1, n
                
        if t1 <= n - sep_b - 2:
            for t2 in range(t1 + sep_b, n-1):
                if y[t2] > y[t1]:
                    
                    dyls=sign_y[t1+1:t2+1,:]
                    xs=x_mat[t1+1:t2+1,:]
                    ixsxs = np.linalg.inv(xs.T @ xs)
                    z = xs.T @ dyls
                    ssr12=(dyls.T @ dyls - z.T @ ixsxs @ z).item()  
                    s2 = (ssr01 + ssr12 + ssrjT_U[t2,0]) / n
                    
                    if s2 < s2_min_UBU:
                        s2_min_UBU, tb_hat_UBU[0], tb_hat_UBU[1] = s2, t1+1, t2+1

                    if y[-1] < y[t2] and t2 < n - 2:
                        s2 = (ssr01 + ssr12 + ssrjT_C[t2,0]) / n
                        if s2 < s2_min_UBC: 
                            s2_min_UBC, tb_hat_UBC[0], tb_hat_UBC[1] = s2, t1+1, t2+1

                    if t1 <= n - sep_b - sep_c - 2 and t2 <= n - sep_c - 2:
                        for t3 in range(t2 + sep_c, n-1):
                            if y[t3] < y[t2]:
                                dyls=sign_y[t2+1:t3+1,:]
                                xs=x_mat[t2+1:t3+1,:]
                                ixsxs = np.linalg.inv(xs.T @ xs)
                                z = xs.T @ dyls
                                ssr23=(dyls.T @ dyls - z.T @ ixsxs @ z).item()
                                
                                s2 = (ssr01 + ssr12 + ssr23 + ssrjT_U[t3,0]) / n
                                if s2 < s2_min_UBCU: 
                                    s2_min_UBCU, tb_hat_UBCU[0], tb_hat_UBCU[1], tb_hat_UBCU[2] = s2, t1+1, t2+1, t3+1

    BIC = np.array([
        n*np.log(s2_min_UB)+2*np.log(n), 
        n*np.log(s2_min_UBU)+3*np.log(n), 
        n*np.log(s2_min_UBC)+4*np.log(n), 
        n*np.log(s2_min_UBCU)+5*np.log(n)])
    
    candidates = np.vstack([tb_hat_UB, tb_hat_UBU, tb_hat_UBC,tb_hat_UBCU])   
    tb_hat_opt = candidates[np.argmin(BIC), :]
    # np.savetxt("series.txt", y, fmt="%.3f")
    return int(tb_hat_opt[0])


def sPSY_df(y):
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    ly = pd.Series(y.flatten()).shift(1).to_numpy().reshape(-1, 1)
    dy = y - ly
    dep = dy[1:, :]
    xmat = ly[1:, :]
    xx = xmat.T @ xmat
    if xx.item() == 0:
        return np.nan
    xxinv = np.linalg.inv(xx)
    b = xxinv @ (xmat.T @ dep)
    r = dep - xmat @ b
    dof = len(xmat) - xmat.shape[1]
    if dof <= 0:
        return np.nan
    sig2 = (r.T @ r).item() / dof
    if sig2 <= 0:
        return np.nan
    df = b / np.sqrt(sig2 * xxinv)
    return df.item()

def sPSY(y):
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    n = len(y)
    r0 = 0.01 + 1.8 / np.sqrt(n)
    w0 = int(np.floor(r0 * n))
    ly = np.vstack([np.array([[np.nan]]), y[:-1, :]])
    dy = y - ly
    y1 = np.divide(dy, np.abs(dy), out=np.zeros_like(dy), where=np.abs(dy) > 0)
    c = np.zeros((n, 1))
    c[1:, :] = np.cumsum(y1[1:, :], axis=0)
    stats = []
    for nr1 in range(n - w0 + 1):
        nr2 = nr1 + w0 - 1
        while nr2 < n:
            stat = sPSY_df(c[nr1:nr2 + 1, :])
            if np.isfinite(stat):
                stats.append(stat)
            nr2 += 1
    return np.max(stats) if stats else np.nan

def weight_one(y, wds):
    y = np.asarray(y, dtype=float).reshape(-1)
    n = len(y)
    x_vec = np.arange(n) / wds
    k_vec = (np.exp(-(x_vec ** 2) / 2) / np.sqrt(2 * np.pi)) * (x_vec > 0) * (x_vec < 1)
    K_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1): K_mat[i, j] = k_vec[i - j]
    row_sums = K_mat.sum(axis=1, keepdims=True)
    return np.divide(K_mat, row_sums, out=np.zeros_like(K_mat), where=row_sums != 0)


def recursive_residual(y, x, d):

    y = np.asarray(y, dtype=float).reshape(-1, 1)
    x = np.asarray(x, dtype=float).reshape(-1, 1)

    n = len(x)

    x_f = np.full((n, 1), np.nan)

    if d == 0:
        x_f = x.copy()

    elif d > 0:
        x_f[d:, :] = x[:-d, :]

    else:
        k = abs(d)
        x_f[:-k, :] = x[k:, :]

    if d >= 0:
        xmat = x_f[d:, :]
        ymat = y[d:, :]
    else:
        xmat = x_f[:-abs(d), :]
        ymat = y[:-abs(d), :]

    xmat = np.column_stack([
        np.ones((len(xmat), 1)),
        xmat
    ])

    k = xmat.shape[1]

    XX = xmat.T @ xmat
    Xy = xmat.T @ ymat

    try:
        beta_temp = np.linalg.solve(XX, Xy)
    except np.linalg.LinAlgError:
        beta_temp = np.linalg.pinv(XX) @ Xy

    res_temp = ymat - xmat @ beta_temp

    return ymat, xmat, res_temp[k:, :]


def recursive_residual2(ymat, xmat, d, f, wds):

    ymat = np.asarray(ymat, dtype=float).reshape(-1, 1)
    xmat = np.asarray(xmat, dtype=float)



    n = len(ymat)
    k1 = xmat.shape[1]

    start_idx = int(f + 1 - abs(d)) - 1
    end_idx = n - 1

    out_len = end_idx - start_idx + 1

    res = np.zeros((out_len, 1))
    sigma_diag = np.zeros((out_len, 1))
    CO_path = np.zeros((out_len, 1))

    CO_sum = 0.0

    r_idx = start_idx
    j = 0

    # Initial OLS
    X_prev = xmat[:r_idx, :]
    y_prev = ymat[:r_idx, :]

    beta_ols_prev = np.linalg.inv(X_prev.T @ X_prev) @ (X_prev.T @ y_prev)

    res_ols_prev = y_prev - X_prev @ beta_ols_prev
    res_ols_prev = res_ols_prev[k1:, :]

    nr_prev = len(res_ols_prev)

    W_prev = weight_one(res_ols_prev, wds)

    sigma_prev = np.sqrt(
        W_prev @ (res_ols_prev ** 2)
    )

    while r_idx <= end_idx:

        left = k1 + wds

        X_slice = xmat[left:r_idx, :]
        y_slice = ymat[left:r_idx, :]

        s_slice = sigma_prev[wds:nr_prev, :]

        usable = min(
            len(X_slice),
            len(y_slice),
            len(s_slice)
        )

        X_slice = X_slice[:usable, :]
        y_slice = y_slice[:usable, :]
        s_slice = s_slice[:usable, :]

        valid = (
            np.isfinite(s_slice[:, 0])
            & (s_slice[:, 0] != 0)
        )

        X0 = X_slice[valid, :] / s_slice[valid, :]
        y0 = y_slice[valid, :] / s_slice[valid, :]

        beta0 = np.linalg.inv(X0.T @ X0) @ (X0.T @ y0) 

        res[j, 0] = (
            ymat[r_idx, 0]
            - (xmat[r_idx:r_idx+1, :] @ beta0).item()
        )

        # Updated OLS
        X_curr = xmat[:r_idx + 1, :]
        y_curr = ymat[:r_idx + 1, :]

        beta_ols_curr =  np.linalg.inv(X_curr.T @ X_curr) @ (X_curr.T @ y_curr) 

        res_ols_curr = y_curr - X_curr @ beta_ols_curr
        res_ols_curr = res_ols_curr[k1:, :]

        nr_curr = len(res_ols_curr)

        W_curr = weight_one(res_ols_curr, wds)

        sigma_curr = np.sqrt(
            W_curr @ (res_ols_curr ** 2)
        )

        sigma_diag[j, 0] = sigma_curr[-1, 0]

        if (
            sigma_diag[j, 0] != 0
            and np.isfinite(sigma_diag[j, 0])
        ):
            CO_sum += res[j, 0] / sigma_diag[j, 0]

        CO_path[j, 0] = abs(CO_sum)

        sigma_prev = sigma_curr
        nr_prev = nr_curr

        j += 1
        r_idx += 1

    return res, CO_path











data = pd.read_csv("index_vs_msft.csv")

t = pd.to_datetime(data.iloc[0:, 0], format="%Y-%m-%d")
# t = (t.astype("int64") // 10**6).to_numpy()

y = pd.to_numeric(data.iloc[0:, 5], errors="coerce").to_numpy()
x = pd.to_numeric(data.iloc[0:, 8], errors="coerce").to_numpy()


ind_d = (np.isfinite(x) & np.isfinite(y) & (x != 0) & (y != 0))

t = t[ind_d]
x = x[ind_d]
y = y[ind_d]


# x = x / np.min(x)
# y = y / np.min(y)
x=x/x[0]
y=y/y[0]

n=len(y)
h = 0.25 * n ** (-0.2)
wds = int(np.floor(h * n))
m=int(0.1*n)
delv = np.arange(-20, 21)






bd = sBIC(x, int(np.floor(0.1*n)), int(np.floor(0.1*n)))
f = int(bd + m)
sPSY_y, sPSY_x = sPSY(y[0:f+1]), sPSY(x[0:f+1])


if sPSY_y > 3.213 and sPSY_x > 3.213:
    CO1 = np.zeros(n)
    CO_temp, ss, ss2 = np.zeros(len(delv)), np.zeros(len(delv)), np.zeros(len(delv))
    
    for k, d0 in enumerate(delv):
        ymat0, xmat0, res_temp0 = recursive_residual(y[:f], x[:f], int(d0))
        W0 = weight_one(res_temp0, wds)
        sigma_hat0 = np.sqrt(W0 @ (res_temp0 ** 2))
        num, den = res_temp0[wds:], sigma_hat0[wds:]
        valid = np.isfinite(den) & (den != 0)
        res0 = num[valid] / den[valid]
        ss[k] = np.mean(res0 ** 2)
        ss2[k] = np.mean(res_temp0[wds:] ** 2)
        CO_temp[k] = abs(np.sum(res0))
    
    ss1_ind, ss2_ind, CO_ind = np.argmin(ss), np.argmin(ss2), np.argmin(CO_temp)
    CO = CO_temp[ss2_ind]
    d = int(delv[ss2_ind])
    cv0 = np.sqrt(4.6) * np.sqrt(f)
    
    if CO < cv0:
        ymat, xmat, res_temp = recursive_residual(y, x, d)
        res, CO_path = recursive_residual2(ymat, xmat, d, f, wds)
        end = min(n, f + len(CO_path))
        CO1[f:end] = CO_path[:end-f,0]

        cv = np.zeros(n)

        for c_g in range(f + 1, n + 1):
            cv[c_g - 1] = np.sqrt(6 + np.log(c_g/f)) * np.sqrt(c_g)

        idx = np.flatnonzero(CO1 > cv)

        if len(idx) == 0: 
            Rejection_data = np.nan
        else:
            Rejection_data = idx[0]
            
        x=pd.Series(x).shift(d)
        
        plt.figure(figsize=(12, 5))
        plt.plot(t, x, label=rf"Nasdaq Computer Index ($\hat{{i}}={d}$)", color="black")
        plt.plot(t, y, label="Microsoft", color="blue")
        plt.xlabel("Time")
        plt.xticks(rotation=45)
        plt.ylabel("Value")
        plt.axvline(x=t.iloc[bd], linestyle="--", label=rf"$\lfloor \hat{{\tau}}_{{1}} T \rfloor$ ({t.iloc[int(bd)].date()})", color="lightgreen")
        plt.axvline(x=t.iloc[f+1], linestyle="--", label=rf"$\lfloor fT \rfloor$ ({t.iloc[int(f+1)].date()})", color="green")
        plt.axvline(x=t.iloc[Rejection_data], linestyle="--", label=rf"Rejection Date ({t.iloc[int(Rejection_data)].date()})", color="red")
        plt.legend()
        plt.tight_layout()
        plt.savefig("coexplosivity_application.png", dpi=180, bbox_inches="tight")
        plt.show()
        print(t.iloc[Rejection_data])
        
    else:
        print("There is no co-bubbling")
        
else:
    print("In either x or y or both x and y, there is no evidence of bubble")
        
  
