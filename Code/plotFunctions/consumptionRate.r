######################################################
## Calculate consumption rate
## using search, attack prob and handling functions
## As well as density and masses
######################################################

V0 = 0.33
D0 = 1.62
H0 = 1
B = 0.75
Rp = 0.1
a = 1
c = 1

calculateSearchRate = function(mi, mj) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*(mj^(0.21))

    return(aij)

}

calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^0.2

    return(Aij)

}


calculateHandling = function(mi, mj) {

    Hij = H0*(mi^(-B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))

    return(Hij)

}

consumptionRate = function(mi, mj, n) {
    
    c = (calculateSearchRate(mi)*calculateAttackProb(mi, mj)*n*mj)/(1 + (calculateSearchRate(mi)*calculateAttackProb(mi, mj)*calculateHandling(mi,mj)*n))

    return(c)

}

dat = cellPopSpec %>% filter(g == max(cellPopSpec$g))
dat1 = expand.grid(1.96, dat$M) %>% as_tibble() %>% rename(mi = 1, M = 2) %>% 
    left_join(dat)
dat1 = dat1 %>% mutate(c = consumptionRate(mi, M, n)*n)
dat2 = dat %>% mutate(c = -consumptionRate(M, 1.96, 1)*n)

sum(dat1$c) + sum(dat2$c)

arrange(dat1, -c)
arrange(dat2, c)

dat3 = expand.grid(1.3, dat$M) %>% as_tibble() %>% rename(mi = 1, M = 2) %>% 
    left_join(dat)
dat3 = dat3 %>% mutate(c = consumptionRate(mi, M, n)*n)
dat4 = dat %>% mutate(c = -consumptionRate(M, 1.3, 2)*n)

sum(dat3$c) + sum(dat4$c)

arrange(dat3, -c)
arrange(dat4, c)
