#library(mco)

#res <- nsga2(hanne4, 2, 2,
#             lower.bounds=c(0, 0), upper.bounds=c(10, 10),
#             constraints=hanne4.constr, cdim=1, generations = 1:100)

#fn <- function(x) c(hanne4(x), hanne4.constr(x))
#objs <- 1:2
#cons <- 3
#ref.point <- rep(10, 2)

#save_run(res, fn = fn, objs = objs, cons = cons, ref.point = ref.point)
