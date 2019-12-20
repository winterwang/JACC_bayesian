#何かのリストを作成
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#stanを実行
fit <- stan(file = '8schools.stan', data = schools_dat)

#少し時間がかかる...
#プロットしてみる
plot(fit)