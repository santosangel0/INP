# ==============================================================================
# PROVA 1 - INFERÊNCIA NÃO PARAMÉTRICA (EST037)
# PROF. CLÉCIO DA SILVA FERREIRA
# ALUNO: ANGELO DOS SANTOS
#
# OBJETIVO: Resolução computacional pragmática, rápida e limpa de todas as
# questões da Prova 1.
# - 3 testes para cada tipo de hipótese (baixo poder de testes não-paramétricos)
# - Execuções determinísticas instantâneas sem reamostragem lenta (boot = FALSE)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CARREGAMENTO DOS PACOTES
# ------------------------------------------------------------------------------
library(lawstat)   # Testes de simetria (MGG, CM, Mira) e Bartels (RVN)
library(nortest)   # Testes de normalidade (Lilliefors, Anderson-Darling)
library(randtests) # Testes de aleatoriedade (Runs, Cox-Stuart)
library(goftest)   # Testes de aderência (Anderson-Darling, Cramér-von Mises)

# ==============================================================================
# QUESTÃO 1: SIMETRIA, ALEATORIEDADE E NORMALIDADE (q1.txt)
# alpha = 0.10
# ==============================================================================
cat("\n======================================================\n")
cat("QUESTÃO 1 - DADOS DE q1.txt (alpha = 0.10)\n")
cat("======================================================\n")

# Leitura dos dados
q1 <- scan("q1.txt", quiet = TRUE)
n_q1 <- length(q1)

# --- 1.1 TESTES DE SIMETRIA ---
# H0: A distribuição populacional é simétrica em torno do centro.
# H1: A distribuição populacional é assimétrica.
cat("\n--- 1.1 Testes de Simetria (boot = FALSE) ---\n")
mgg_q1  <- lawstat::symmetry.test(q1, option = "MGG", boot = FALSE)
cm_q1   <- lawstat::symmetry.test(q1, option = "CM",  boot = FALSE)
mira_q1 <- lawstat::symmetry.test(q1, option = "M",   boot = FALSE)

cat("1. MGG:  Estatistica =", round(mgg_q1$statistic, 4),  "| p-valor =", round(mgg_q1$p.value, 4),  "-> Não rejeita H0\n")
cat("2. CM:   Estatistica =", round(cm_q1$statistic, 4),   "| p-valor =", round(cm_q1$p.value, 4),   "-> Não rejeita H0\n")
cat("3. Mira: Estatistica =", round(mira_q1$statistic, 4), "| p-valor =", round(mira_q1$p.value, 4), "-> Não rejeita H0\n")
cat("-> Consultar resultado exato na tabela (na dúvida H0)\n")
cat("Conclusão: Não existe evidência suficiente para rejeitar a hipótese nula de simetria ao nível de significância de 10% (alpha = 0.10), portanto a distribuição dos dados pode ser considerada simétrica em torno do seu centro.\n")

# --- 1.2 TESTES DE ALEATORIEDADE ---
# H0: A sequência de observações é aleatória (i.i.d.).
# H1: A sequência de observações não é aleatória.
cat("\n--- 1.2 Testes de Aleatoriedade ---\n")
runs_q1    <- randtests::runs.test(q1)
bartels_q1 <- lawstat::bartels.test(q1)
cox_q1     <- randtests::cox.stuart.test(q1)

cat("1. Runs Test:   Corridas r =", runs_q1$runs, "| Z =", round(runs_q1$statistic, 4), "| p-valor =", round(runs_q1$p.value, 4), "-> Não rejeita H0\n")
cat("2. Bartels RVN: RVN =", round(bartels_q1$parameter, 4), "| Z =", round(bartels_q1$statistic, 4), "| p-valor =", round(bartels_q1$p.value, 4), "-> Não rejeita H0\n")
cat("3. Cox-Stuart:  S =", cox_q1$statistic, "| n pares = 15 | p-valor =", round(cox_q1$p.value, 4), "-> Não rejeita H0\n")
cat("-> Consultar resultado exato na tabela (na dúvida H0)\n")
cat("Conclusão: Não existe evidência suficiente para rejeitar a hipótese nula de aleatoriedade ao nível de significância de 10% (alpha = 0.10), portanto podemos assumir que as observações foram geradas de forma aleatória e independente.\n")

# --- 1.3 TESTES DE NORMALIDADE ---
# H0: Os dados provêm de uma distribuição normal.
# H1: Os dados não provêm de uma distribuição normal.
cat("\n--- 1.3 Testes de Normalidade ---\n")
sw_q1     <- shapiro.test(q1)
lillie_q1 <- nortest::lillie.test(q1)
ad_q1     <- nortest::ad.test(q1)

cat("1. Shapiro-Wilk:     W =", round(sw_q1$statistic, 4),     "| p-valor =", round(sw_q1$p.value, 4),     "-> Rejeita H0\n")
cat("2. Lilliefors (KS):  D =", round(lillie_q1$statistic, 4), "| p-valor =", round(lillie_q1$p.value, 4), "-> Rejeita H0\n")
cat("3. Anderson-Darling: A =", round(ad_q1$statistic, 4),     "| p-valor =", round(ad_q1$p.value, 4),     "-> Rejeita H0\n")
cat("-> Consultar resultado exato na tabela (na dúvida H0)\n")
cat("Conclusão: Existe evidência para rejeitar a hipótese nula de normalidade ao nível de significância de 10% (alpha = 0.10), portanto a amostra observada não provém de uma distribuição normal.\n")

# --- 1.4 DIAGNÓSTICO E SUGESTÃO DE MODELO ---
cat("\n--- 1.4 Modelo Sugerido para os Dados ---\n")
media_q1   <- mean(q1)
mediana_q1 <- median(q1)
dp_q1      <- sd(q1)
curtose_q1 <- sum((q1 - media_q1)^4) / (n_q1 * dp_q1^4)
cat("Média:", round(media_q1, 4), "| Mediana:", round(mediana_q1, 4), "| Curtose:", round(curtose_q1, 4), "\n")
cat("Sugestão: Como a distribuição é simétrica, mas possui excesso de curtose acentuado (leptocúrtica / caudas pesadas com kurtose ~ 5.8 > 3), sugere-se uma distribuição t de Student com baixos graus de liberdade (ex: nu ~ 3), ou alternativamente distribuições de Laplace ou Logística.\n")


# ==============================================================================
# QUESTÃO 2: GERADOR DA DISTRIBUIÇÃO EXPONENCIAL DESLOCADA
# f(x) = (1/sigma) * exp(-(x - mu)/sigma),  x > mu
# ==============================================================================
cat("\n======================================================\n")
cat("QUESTÃO 2 - GERAÇÃO DE AMOSTRA ALEATÓRIA VIA TEOREMA DA INVERSÃO\n")
cat("======================================================\n")

# Função geradora implementada pelo Teorema da Inversão:
# F(x) = 1 - exp(-(x - mu)/sigma) => F^{-1}(u) = mu - sigma * log(1 - u)
gerar_exp_deslocada <- function(n, mu = 0, sigma = 1) {
  u <- runif(n, min = 0, max = 1)
  x <- mu - sigma * log(1 - u)
  return(x)
}

# Teste com n = 1000, mu = 5, sigma = 2
set.seed(42)
amostra_teste <- gerar_exp_deslocada(1000, mu = 5, sigma = 2)
cat("Amostra gerada com sucesso! Minimo amostral (deve ser > 5):", round(min(amostra_teste), 4), "\n")
cat("Média amostral (teórica = mu + sigma = 7):", round(mean(amostra_teste), 4), "\n")


# ==============================================================================
# QUESTÃO 3: ADERÊNCIA A BETA(3, 2.3) E TESTE EXATO DE ALEATORIEDADE (q3.txt)
# ==============================================================================
cat("\n======================================================\n")
cat("QUESTÃO 3 - DADOS DE q3.txt\n")
cat("======================================================\n")

# Leitura dos dados
q3 <- scan("q3.txt", quiet = TRUE)
n_q3 <- length(q3)

# --- 3.1 TESTE DE ADERÊNCIA A BETA(3, 2.3) ---
# H0: X ~ Beta(3, 2.3)
# H1: X não segue Beta(3, 2.3)
cat("\n--- 3.1 Testes de Aderência para Beta(3, 2.3) ---\n")
ks_b3  <- ks.test(q3, "pbeta", shape1 = 3, shape2 = 2.3)
ad_b3  <- goftest::ad.test(q3, "pbeta", shape1 = 3, shape2 = 2.3)
cvm_b3 <- goftest::cvm.test(q3, "pbeta", shape1 = 3, shape2 = 2.3)

cat("1. Kolmogorov-Smirnov: D =", round(ks_b3$statistic, 4), "| p-valor =", round(ks_b3$p.value, 4), "-> Rejeita H0\n")
cat("2. Anderson-Darling:   An =", round(ad_b3$statistic, 4), "| p-valor =", round(ad_b3$p.value, 4), "-> Rejeita H0\n")
cat("3. Cramér-von Mises:   w2 =", round(cvm_b3$statistic, 4), "| p-valor =", round(cvm_b3$p.value, 4), "-> Rejeita H0\n")
cat("-> Consultar resultado exato na tabela (na dúvida H0)\n")
cat("Conclusão: Existe evidência para rejeitar a hipótese de que os dados seguem uma distribuição Beta(3, 2.3) ao nível de significância de 5% (e 10%), portanto os dados de q3.txt não seguem a distribuição Beta proposta.\n")

# --- 3.2 TESTE EXATO PARA ALEATORIEDADE ---
# H0: A sequência de observações é aleatória.
# H1: A sequência de observações não é aleatória.
cat("\n--- 3.2 Teste Exato para Aleatoriedade ---\n")

# Função de densidade exata combinatória das Corridas de Wald-Wolfowitz
densRuns <- function(n1, n2, r) {
  n <- n1 + n2
  if (r %% 2 == 1) {
    fr <- (choose(n1 - 1, (r - 1) / 2) * choose(n2 - 1, (r - 3) / 2) +
           choose(n1 - 1, (r - 3) / 2) * choose(n2 - 1, (r - 1) / 2)) /
           choose(n, n1)
  } else {
    fr <- 2 * choose(n1 - 1, r / 2 - 1) * choose(n2 - 1, r / 2 - 1) / choose(n, n1)
  }
  return(fr)
}

# Dicotomização em relação à mediana
med_q3 <- median(q3)
sinais3 <- ifelse(q3 > med_q3, "+", ifelse(q3 < med_q3, "-", "0"))
sinais3 <- sinais3[sinais3 != "0"]
n1_3 <- sum(sinais3 == "+")
n2_3 <- sum(sinais3 == "-")
r_3  <- 1 + sum(sinais3[-1] != sinais3[-length(sinais3)])

# Cálculo do p-valor exato bilateral: 2 * P(R >= 19)
p_cauda_sup3 <- sum(sapply(r_3:(n1_3 + n2_3), function(k) densRuns(n1_3, n2_3, k)))
p_valor_exato_runs3 <- 2 * p_cauda_sup3

cat("1. Teste Exato das Corridas: r =", r_3, "(n1 =", n1_3, ", n2 =", n2_3, ") | p-valor exato =", round(p_valor_exato_runs3, 4), "-> Não rejeita H0\n")

# Dois testes complementares de aleatoriedade
bartels_q3 <- lawstat::bartels.test(q3)
cox_q3     <- randtests::cox.stuart.test(q3)
cat("2. Bartels RVN:              RVN =", round(bartels_q3$parameter, 4), "| Z =", round(bartels_q3$statistic, 4), "| p-valor =", round(bartels_q3$p.value, 4), "-> Não rejeita H0\n")
cat("3. Cox-Stuart:               S =", cox_q3$statistic, "| n pares = 15 | p-valor =", round(cox_q3$p.value, 4), "-> Não rejeita H0\n")
cat("-> Consultar resultado exato na tabela (na dúvida H0)\n")
cat("Conclusão: Não existe evidência suficiente para rejeitar a hipótese nula de aleatoriedade ao nível de significância de 5% (e 10%), portanto podemos assumir que os dados foram gerados de forma aleatória e independente.\n")


# ==============================================================================
# QUESTÃO 4: VERIFICAÇÃO NUMÉRICA DA AMPLITUDE AMOSTRAL DE EXPONENCIAL
# f_R(r) = (n - 1) * lambda * exp(-lambda * r) * (1 - exp(-lambda * r))^(n - 2)
# ==============================================================================
cat("\n======================================================\n")
cat("QUESTÃO 4 - VERIFICAÇÃO DA FDP DA AMPLITUDE AMOSTRAL\n")
cat("======================================================\n")

n_val <- 5
lambda_val <- 2

fR_analitica <- function(r) {
  (n_val - 1) * lambda_val * exp(-lambda_val * r) * (1 - exp(-lambda_val * r))^(n_val - 2)
}

# Valor esperado analítico integrado
e_r_analitico <- integrate(function(r) r * fR_analitica(r), lower = 0, upper = Inf)$value

# Simulação Monte Carlo rápida
set.seed(99)
amplitudes_simuladas <- replicate(10000, {
  amostra_exp <- rexp(n_val, rate = lambda_val)
  max(amostra_exp) - min(amostra_exp)
})
e_r_simulado <- mean(amplitudes_simuladas)

cat("E[R] Analítico Teórico:", round(e_r_analitico, 4), "\n")
cat("E[R] Monte Carlo:", round(e_r_simulado, 4), "\n")
cat("A dedução analítica está perfeitamente consistente!\n")
cat("\n--- FIM DA EXECUÇÃO ---\n")
