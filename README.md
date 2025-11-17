##### dados #####

dados_env <- data.frame(
  site = c("T1","T2","T3","T4","T5","T9","R1"),
  riqu_core = c(5,0,2,3,3,0,2),
  abund_core = c(16,0,5,6,37,0,2),
  riqu_sapo = c(4,8,7,4,4,2,8),
  abund_sapo = c(18,26,21,19,8,8,28),
  freq_bromelia = c(1,0.52,0.56,1,0.12,0.76,1),
  qtd_rio = c(0,2,1,0,0,3,3),
  dist_rio = c(284.049,92.436,121.875,168.458,175.801,90.356,0),
  umid = c(80.95,84,97.5,85.2,81.95,84.15,80.95),
  temp = as.numeric(c(26.2,22.4,23.15,23.05,23.8,23.3,25.05))
)


############ 1ª etapa - Análises exploratórias #############

library(MASS)

# riqu_core <- dados_env$riqu_core
# abund_core <- dados_env$abund_core
# bromelia <- dados_env$freq_bromelia

##### PONCHO ####

#Function

poncho<-function(
    table,col=1,col.places=2,col.species=3,col.gradient="gray",
    highlight.species=NA,highlight.places=NA,
    xyratio=1.3333333333,cex.species=NULL,lwd=1,
    cex.lab=NULL,cex.title=3,cex.gradient=2,gradient=NULL,
    file=NULL,title=NULL,xlab="Ordered places",
    ylab="Abundance",show.grad=TRUE,show.places=NULL,
    xfactor=0,gradlab="Gradient",lty=1,lty.lines=1,border=1,decrease=FALSE,add=F,spgradient=NULL)
{
  
  if(is.null(file)){
    if(is.null(cex.species)){cex.species=.5}
    if(is.null(cex.lab)){cex.lab=2}
    if(is.null(cex.gradient)){cex.gradient=1.5}
  }
  
  if(is.null(cex.gradient)){cex.gradient=2}
  if(is.null(cex.lab)){cex.lab=3}
  
  ####################################
  
  if(is.null(file)==FALSE){
    
    pdf(file=file,,paper="A4",width=0,height=0)
    
  }
  
  par(mar=c(2,3,.8,3))
  
  #################################
  
  if(is.null(gradient)){gradient<-prcomp(table)$x[,1]}
  
  multiplier=ifelse(decrease==F,1,-1)
  
  
  if(is.null(spgradient)==F){
    
    tabela<-table[order(multiplier*gradient),order(spgradient)]
    if(length(highlight.species)==length(col.species)){
      col.species<-col.species[order(spgradient)]
    }
  }else{
    
    tabela<-table[order(multiplier*gradient),order(colSums(table*multiplier*gradient)/colSums(table))]
    
    if(length(highlight.species)==length(col.species)){
      col.species<-col.species[order(colSums(table*multiplier*gradient)/colSums(table))]
    }
    
    
  }
  
  
  #################################
  x<-90
  y<-x*xyratio
  
  dimx<-x/length(0:(nrow(tabela)-1))
  dimy<-y/length(1:ncol(tabela))
  showg=0
  if(is.null(show.grad)==F){showg=22}
  
  if(add==F){
    plot.new()
    plot.window(xlim=c(0,120+xfactor),ylim=c(0,120+showg))
  }
  
  col<-rep(col,nrow(tabela))
  col[match(highlight.places,rownames(tabela))]<-col.places
  
  for(i in 1:ncol(tabela)){
    
    col.sp<-match(highlight.species,colnames(tabela)[i])
    
    if(sum(col.sp,na.rm=T)!=0){col.2=col.species}else{col.2=col}
    
    if(length(highlight.species)==length(col.species)){col.2=col.species[i]}
    
    rect((0:(nrow(tabela)-1))*dimx,(i-1)*dimy,(1:nrow(tabela))*dimx-dimx*.1,(i-1)*dimy+(tabela[,i]/max(tabela))*dimy*.8,col=col.2,lty=lty,lwd=lwd,border=border)
    lines(c(-1,(nrow(tabela))*dimx),c((i-1)*dimy,(i-1)*dimy),lty=lty.lines,lwd=lwd)
  }
  text(rep(x+2,ncol(tabela)),cumsum(rep(120/ncol(tabela),ncol(tabela)))-120/ncol(tabela)/2,gsub("_"," ",colnames(tabela)),cex=if(is.null(cex.species)){ifelse(30/ncol(tabela)>1,1,30/ncol(tabela))}else{cex.species},adj=0,font=3)
  
  if(lty.lines!=0){
    rect(-1,0,-1,(ncol(tabela)-1)*dimy+dimy*.8,lwd=lwd)
  }
  
  text(x+4,133.5,gradlab,cex=cex.gradient,adj=0)
  
  par(mgp=c(1,4,5))
  title(main=title,cex.main=cex.title)
  
  par(mgp=c(0,4,5))
  title(ylab=ylab,cex.lab=cex.lab)
  
  par(mgp=c(1,8,10),mar=c(3.5,3,.8,8))
  title(xlab=xlab,cex.lab=cex.lab)
  
  
  par(mgp=c(4,8,10),mar=c(0,0,0,0))
  if(is.null(show.grad)==F){
    gradient<-sort(gradient,decreasing=decrease)
    gradient2<-gradient/abs(max(c(gradient,0))-min(c(gradient,0)))
    val=126-min(c(gradient2*15,0))
    
    rect((0:(nrow(tabela)-1))*dimx,val,(1:nrow(tabela))*dimx-dimx*.1,val+gradient2*15,lwd=lwd,col=col.gradient)
    
    rect(x+2,126,x+2,126+15,lwd=lwd)
    
    rect(x+2,126,x+3,126,lwd=lwd)
    
    rect(x+2,126+15,x+3,126+15,lwd=lwd)
    
    text(x+3.5,126,round(min(c(0,gradient)),2),adj=0,cex=.7)
    text(x+3.5,126+15,round(max(c(0,gradient)),2),adj=0,cex=.7)
    
  }
  
  par(mar=c(5, 4, 4, 2) + 0.1,mgp=c(0,-.5,0))
  if(is.null(show.places)==F){
    axis(1,at=((0:(nrow(tabela)-1))*dimx+(1:nrow(tabela))*dimx-dimx*.1)/2,labels=rownames(tabela),las=2,cex.axis=if(is.null(cex.species)){ifelse(30/ncol(tabela)>1,1,30/ncol(tabela))}else{cex.species},tick=F)
  }
  #title(main="B",adj=0,cex.main=2)
  par(mar=c(5, 4, 4, 2) + 0.1,mgp=c(3, 1, 0))
  if(is.null(file)==F){
    
    dev.off()
    paste("You have saved the image as",file)
  }else{"If resolution is low, try to use the argument file=imagename.pdf to save the image in a file"}
  
}


#end

##Use
#poncho(table)

install.packages("tidyverse")  
library(tidyverse) 

########### AQUI COMEÇA ####################


# Criando um data frame com os dados de abundancia de anfíbios (pode fazer o mesmo para core) 
data <- data.frame(
  site = c("T1", "T2", "T3", "T4", "T5", "T9", "R1"),
  Adenomera_araucaria = c(0, 0, 1, 0, 0, 0, 1),
  Aplastodiscus_perviridis = c(0, 0, 0, 0, 0, 0, 3),
  Boana_bischoffi = c(0, 0, 0, 0, 1, 0, 0),
  Boana_faber = c(2, 0, 0, 0, 0, 0, 0),
  Bokermannohyla_hylax = c(0, 0, 2, 2, 4, 0, 8),
  Dendroprhyniscus_berthalutzae = c(1, 5, 7, 5, 2, 1, 2),
  Fritziana_mitus = c(10, 11, 6, 11, 0, 0, 8),
  Hylodes_meridionalis = c(0, 0, 0, 0, 0, 0, 2),
  Ischnocnema_henseli = c(5, 5, 2, 0, 0, 7, 3),
  Ischnocnema_aff._manezinho = c(0, 1, 2, 0, 0, 0, 0),
  Phyllomedusa_distincta = c(0, 1, 0, 0, 0, 0, 0),
  Physalaemus_lateristriga = c(0, 0, 1, 0, 0, 0, 0),
  Physalaemus_nanus = c(0, 0, 0, 0, 0, 0, 1),
  Proceratophrys_boiei = c(0, 1, 0, 0, 1, 0, 0),
  Rhinella_ornata = c(0, 1, 0, 0, 0, 0, 0),
  Rhinella_icterica = c(0, 0, 0, 1, 0, 0, 0),
  Scinax_catharinae = c(0, 1, 0, 0, 0, 0, 0),
  dist_rio = c(284.049, 92.436, 121.875, 168.458, 175.801, 90.356, 0)
)


# convertendo dados para o formato "wide para ggplot
data_long <- data %>%
  pivot_longer(cols = -c(site, dist_rio), 
               names_to = "species", 
               values_to = "abundance") %>%
  filter(abundance > 0)  # Remove species with zero abundance

# Stacked bar plot
ggplot(data_long, aes(x = reorder(site, dist_rio), y = abundance, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Amphibian Species Abundance per Site (Ordered by River Distance)",
    x = "Site (Ordered by Distance from River)",
    y = "Abundance",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(labels = function(x) gsub("_", " ", x))

########################## PONCHO ###############################

############## ANUROS ###########

# dados abundancia
data_matrix <- matrix(c(
  0, 0, 0, 2, 0, 1, 10, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 5, 11, 0, 5, 1, 1, 0, 0, 1, 1, 0, 1,
  1, 0, 0, 0, 2, 7, 6, 0, 2, 2, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 2, 5, 11, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 1, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 3, 0, 0, 8, 2, 8, 2, 3, 0, 0, 0, 1, 0, 0, 0, 0
), nrow = 7, byrow = TRUE)

# nome das linhas e colunas
rownames(data_matrix) <- c("T1", "T2", "T3", "T4", "T5", "T9", "R1")
colnames(data_matrix) <- c("A._araucaria", "A._perviridis", "B._bischoffi",
                           "B._faber", "B._hylax", "D._berthalutzae",
                           "F._mitus", "H._meridionalis", "I._henseli",
                           "I._aff._manezinho", "P._distincta", 
                           "P._lateristriga", "P._nanus", "P._boiei",
                           "R._ornata", "R._icterica", "S._catharinae")

# gradiente de distância (pode substituir por qualquer outra variável contínua)
river_gradient <- c(284.049, 92.436, 121.875, 168.458, 175.801, 90.356, 0)
freq_bromelia <- c(1.0, 1.0, 1,0, 0.76, 0,56, 0,52)


# função básica
poncho(table = data_matrix, 
       gradient = river_gradient,
       xlab = "Pontos",
       ylab = "Abundância",
       gradlab = "Distância do Rio (m)",
       show.places = TRUE,
       cex.species = 3.0,
       cex.lab = 3.5,
       col.gradient = "lightblue")



# Ordenando por espécies-alvo usando o parâmetro spgradient
# Por exemplo, ordenar pela abundância total:

species_order <- order(colSums(data_matrix))

poncho(table = data_matrix,
       gradient = river_gradient,
       xlab = "Pontos",
       ylab = "Abundância",
       gradlab = "Distância do Rio (m)",
       cex.species = 0.6,
       cex.lab = 1.0,
       spgradient = species_order,
       title = "Amphibian Distribution - Species Ordered by Total Abundance")

################ CORETHRELLA #################

library(vegan)

# Matriz de abundâncias
data_core <- matrix(c(
  1, 6, 4, 0, 3, 2,  # T1
  0, 0, 0, 0, 0, 0,  # T2
  4, 0, 1, 0, 0, 0,  # T3
  1, 1, 4, 0, 0, 0,  # T4
  5,10,22, 0, 0, 0,  # T5
  0, 0, 0, 0, 0, 0,  # T9
  1, 0, 0, 1, 0, 0   # R1
), nrow = 7, byrow = TRUE)

rownames(data_core) <- c("T1","T2","T3","T4","T5","T9","R1")
colnames(data_core) <- c("C._lopesi","C._cardosoi","C._pillosa","C._atricornis","C._blanda","C._alticola")

# Gradiente de distância
river_gradient_core <- c(284.049, 92.436, 121.875, 168.458, 175.801, 90.356, 0)

# Ordenar sites pelo gradiente
order_sites_core <- order(river_gradient_core, decreasing = TRUE)
data_core_ord <- data_core[order_sites_core, ]
river_gradient_core_ord <- river_gradient_core[order_sites_core]

# Plot poncho final com nomes dos sites horizontalmente

poncho (table = data_core_ord,
        gradient = river_gradient_core_ord,
        gradlab = "Distância do Rio (m)",
        xlab = "Pontos",
        ylab = "Abundância",
        cex.species = 0.7,  
        cex.lab = 1.0)
poncho         


species_order <- order(colSums(data_core))

poncho(table = data_core,
       gradient = river_gradient,
       xlab = "Pontos",
       ylab = "Abundância",
       gradlab = "Distância do Rio (m)",
       cex.species = 1.0,
       cex.lab = 1.2,
       spgradient = species_order)

freq_bromelia <- c(1, 0.52, 0.56, 1, 0.12, 0.76, 1)

poncho(table = data_matrix,
       gradient = freq_bromelia,
       xlab = "Pontos",
       ylab = "Abundância",
       gradlab = "Índice de Bromélia",
       cex.species = 0.6,
       col.gradient = "lightblue",
       cex.lab = 0.8)

poncho(table = data_core,
       gradient = freq_bromelia,
       xlab = "Pontos",
       ylab = "Abundância",
       gradlab = "Índice de Bromélia",
       cex.species = 0.6,
       col.gradient = "lightblue",
       cex.lab = 0.8)


nrow(data_matrix)
length(freq_bromelia)

dim(data_core)
length(freq_bromelia)
length(species_order)

######## Gráficos de Dispersão ######

install.packages("ggplot2")
library(ggplot2)

#Gráfico de dispersão riqueza X rio mais próximo
ggplot(dados_env, aes(x = dist_rio, y = riqu_core)) +
  geom_point(size = 4, alpha = 0.7) +        # pontos transparentes e maiores
  geom_smooth(method = "lm") +                # linha de regressão linear ajustada
  labs(x = "Distância do rio mais próximo", 
       y = "Riqueza observada") +             # nomes dos eixos
  theme_minimal()            

ggplot(dados_env, aes(x = freq_bromelia, y = riqu_core)) +
  geom_point(size = 4, alpha = 0.7) +        # pontos transparentes e maiores
  geom_smooth(method = "lm") +                # linha de regressão linear ajustada
  labs(x = "Bromelia", 
       y = "Riqueza observada") +             # nomes dos eixos
  theme_minimal()   

#Gráfico de dispersão abund X rio mais próximo
ggplot(dados_env, aes(x =dist_rio, y = abund_core)) +
  geom_point(size = 4, alpha = 0.7) +        # pontos transparentes e maiores
  geom_smooth(method = "lm") +                # linha de regressão linear ajustada
  labs(x = "Distância do rio mais próximo", 
       y = "Abundância observada") +             # nomes dos eixos
  theme_minimal()            


ggplot(dados_env, aes(x = dist_rio, y = abund_sapo)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(title = "Distancia do rio  vs abundancia sapo")

ggplot(dados, aes(x = freq_bromelia, y = Valor, fill = Composicao)) +
  geom_col(position = "dodge") +         # barras lado a lado para comparação
  labs(x = "Índice de Bromélia", y = "Composição de Corethrella") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# distribuição das variáveis respostas = riq e abund

plot(dados_env$riqu_core)

########### Colinearidade variáveis explicativas = todas ############

# calcular matriz de correlação
matriz_cor <- cor(dados_new, use = "complete.obs")

# plot da correlação
corrplot::corrplot(matriz_cor)

##Resultado: dist_rio e qtd_rio/abund_sapo e riqu_sapo correlacionadas

################################# Modelagem ########################################

# 2ª etapa - Modelo Global - com TODAS as variáveis explicativas que sobraram depois da análise de colinearidade - Um modelo para riqueza e outro para abundância
# 3ª etapa - Drop

##Riqueza 

model_riq1 <- glm(riqu_core ~ abund_sapo + riqu_sapo + freq_bromelia + dist_rio + umid + temp, family = poisson, data = dados_env)

model_riq2 <- glm(riqu_core ~ abund_sapo + freq_bromelia + dist_rio + umid + temp, family = poisson, data = dados_env)

model_riq3 <- glm(riqu_core ~ freq_bromelia + dist_rio + umid + temp, family = poisson, data = dados_env)

model_riq4 <- glm(riqu_core ~ dist_rio + umid + temp, family = poisson, data = dados_env)

model_riq5 <- glm(riqu_core ~ temp + dist_rio, family = poisson, data = dados_env)

model_riq6 <- glm(riqu_core ~ dist_rio, family = poisson, data = dados_env) #ESSE

summary(model_riq6)
plot(model_riq6)


##Abundancia


model_abund1 <- glm(abund_core ~ abund_sapo + freq_bromelia + dist_rio + umid + temp, family = poisson, data = dados_env)

model_abund2 <- glm(abund_core ~ abund_sapo + freq_bromelia + dist_rio + temp, family = poisson, data = dados_env)

model_abund3 <- glm(abund_core ~ freq_bromelia + dist_rio + temp, family = poisson, data = dados_env)

model_abund4 <- glm(abund_core ~ freq_bromelia + dist_rio, family = poisson, data = dados_env) #ESSE

model_abund5 <- glm(abund_core ~ freq_bromelia + dist_rio, family = gaussian(), data = dados_env)

# Ajustar modelo binomial negativo, que estima theta automaticamente
model_nb <- glm.nb(abund_core ~ freq_bromelia + dist_rio, data = dados_env)

summary(model_abund4)
summary(model_riq6)

# 4ª etapa - Seleção de Modelos - AIC

AIC(model_riq1, model_riq2, model_riq3, model_riq4, model_riq5, model_riq6)

AIC(model_abund1, model_abund2, model_abund3, model_abund4, model_abund5, model_nb)

#####Escolhidos: model_abund4 e model_riqu6

############Gráfico com valores estimados para abundância X dist do rio#########

novo_dist <- data.frame(
  dist_rio = seq(min(dados_env$dist_rio), max(dados_env$dist_rio), length.out = 100),
  freq_bromelia = mean(dados_env$freq_bromelia)
)

# previsão com erro padrão
pred_dist <- predict(model_abund4, newdata = novo_dist, type = "link", se.fit = TRUE)

# valores previstos na escala da resposta
novo_dist$fit <- exp(pred_dist$fit)
novo_dist$lwr <- exp(pred_dist$fit - 1.96 * pred_dist$se.fit)
novo_dist$upr <- exp(pred_dist$fit + 1.96 * pred_dist$se.fit)

ggplot(dados_env, aes(x = dist_rio, y = abund_core)) +
  geom_point(alpha = 0.6) +
  geom_line(data = novo_dist, aes(x = dist_rio, y = fit), 
            color = "red", size = 1, inherit.aes = FALSE) +
  geom_ribbon(data = novo_dist, aes(x = dist_rio, ymin = lwr, ymax = upr), 
              inherit.aes = FALSE, alpha = 0.2, fill = "red") +
  theme_minimal() +
  labs(x = "Distância do rio", y = "Abundância de Corethrella")

############Gráfico com valores estimados para abundância X freq de bromelia#####

predito <- predict(model_abund4, type = "response")

ggplot(dados_env, aes(x = freq_bromelia)) +
  geom_point(aes(y = abund_core), size = 4, shape = 21, fill = "gray", alpha = 0.7) +
  geom_line(data = novos_dados, aes(x = freq_bromelia, y = predito), size = 1) +
  labs(x = "Frequência de bromélia", y = "Abundância prevista de Corethrella") +
  theme_minimal()

# criar grid de valores para freq_bromelia
novo_brom <- data.frame(
  freq_bromelia = seq(min(dados_env$freq_bromelia),
                      max(dados_env$freq_bromelia),
                      length.out = 100),
  dist_rio = mean(dados_env$dist_rio) # fixa a outra variável no valor médio
)

# gerar previsões com erro padrão na escala do link (log)
pred_brom <- predict(model_abund4, newdata = novo_brom, type = "link", se.fit = TRUE)

# converter para escala da resposta (contagem esperada)
novo_brom$fit <- exp(pred_brom$fit)
novo_brom$lwr <- exp(pred_brom$fit - 1.96 * pred_brom$se.fit)
novo_brom$upr <- exp(pred_brom$fit + 1.96 * pred_brom$se.fit)

# gráfico
ggplot(dados_env, aes(x = freq_bromelia, y = abund_core)) +
  geom_point(alpha = 0.6) +  # dados observados
  geom_line(data = novo_brom, aes(x = freq_bromelia, y = fit), 
            color = "blue", size = 1, inherit.aes = FALSE) +
  geom_ribbon(data = novo_brom, aes(x = freq_bromelia, ymin = lwr, ymax = upr), 
              inherit.aes = FALSE, alpha = 0.2, fill = "blue") +
  theme_minimal() +
  labs(x = "Frequência de bromélias", 
       y = "Abundância de Corethrella")


############Gráfico com valores estimados para riqueza X dist do rio#########

library(ggplot2)


predito <- predict(model_riq6, type = "response")

ggplot(data = dados_env, aes(x = dist_rio)) +
  geom_point(aes(y = riqu_core), size = 2.5, shape = 21, fill = "gray", alpha = 0.7) +
  geom_line(data = dados_env, aes(x = dist_rio, y = predito), size = 1) +
  labs(x = "Distância do rio", y = "Riqueza prevista de Corethrella") +
  theme_minimal()


library(ggplot2)

# criar grid de valores para dist_rio
novo_riq <- data.frame(
  dist_rio = seq(min(dados_env$dist_rio),
                 max(dados_env$dist_rio),
                 length.out = 100)
)

# gerar previsões com erro padrão na escala do link (log)
pred_riq <- predict(model_riq6, newdata = novo_riq, type = "link", se.fit = TRUE)

# converter para escala da resposta (abundância esperada)
novo_riq$fit <- exp(pred_riq$fit)
novo_riq$lwr <- exp(pred_riq$fit - 1.96 * pred_riq$se.fit)
novo_riq$upr <- exp(pred_riq$fit + 1.96 * pred_riq$se.fit)

# gráfico
ggplot(dados_env, aes(x = dist_rio, y = riqu_core)) +
  geom_point(alpha = 0.6) +  # dados observados
  geom_line(data = novo_riq, aes(x = dist_rio, y = fit), 
            color = "green", size = 1, inherit.aes = FALSE) +
  geom_ribbon(data = novo_riq, aes(x = dist_rio, ymin = lwr, ymax = upr), 
              inherit.aes = FALSE, alpha = 0.2, fill = "green") +
  theme_minimal() +
  labs(x = "Distância do rio", 
       y = "Riqueza de Corethrella")



##### 5ª etapa - Validação de Modelos


#Calculo pseudo R quadrado

install.packages("pscl")
library(pscl)

pR2(model_abund4) #0.68
pR2(model_riq6) #0.17



###### Curva de rarefação ########

install.packages("iNEXT")
library(iNEXT)

head(dados_rarefacao)

dados_rarefacao <- data.frame(
  Ambrozio_2019 = c(0,0,0,0,0,0,0,0,0,0,0,53,25,180,117,3,3,1,0,0,0,0,0,0,0,0,0,0),
  Caldart_2016  = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,68,0,0,0,82,35,2,1,0,0,0,0,0,0),
  Ambrozio_2012 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,15,2,2,2),
  This_work     = c(0,0,0,0,0,0,0,0,0,0,0,0,0,31,12,0,3,0,0,2,1,0,17,0,0,0,0,0),
  row.names = c(
    "C_menini", "C_amazonica", "C_unifasciata", "C_anniae", "C_bifida_davisi",
    "C_appendiculata", "C_quadrivittata", "C_orthicola", "C_edwardsi", "C_brandiae",
    "C_yanomami", "C_inca", "C_carariensis", "C_pillosa", "C_lopesi",
    "C_amabils", "C_blanda", "C_infuscata", "C_yucuman", "C_alticola",
    "C_atricornis", "C_wirthi", "C_cardosoi", "sp_1", "sp_2", "sp_3", "sp_4", "sp_5"
  )
)

## Número de indivíduos por local
colSums(dados_rarefacao)


# Análise 
resultados_rarefacao <- iNEXT(dados_rarefacao, q = 0, datatype = "abundance", endpoint = 800)

# Gráfico
ggiNEXT(resultados_rarefacao, type = 1) +
  geom_vline(xintercept = 75, lty = 2) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  labs(x = "Número de indivíduos", y = "Riqueza de espécies") +
  theme_minimal()

############### ajustado ###################

g_ludico <- ggplot(df_plot, aes(x = x, y = y, color = site)) +
  # Linha de referência
  geom_vline(xintercept = 66, linetype = 2, color = "gray40", linewidth = 0.7, alpha = 0.8) +
  scale_alpha_manual(values = c(0.7, 0.4) ) +  # 0.7 para linhas, 0.4 para áreas (opcional)
  labs(
    title = "Curvas de rarefação da riqueza de Corethrella em diferentes estudos",
    x = "Número de indivíduos",
    y = "Riqueza de espécies",
    linetype = "Tipo de curva"
  ) +
  theme_minimal(base_family = "nunito") +
  theme(
    text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.background = element_rect(fill = "white", color = NA)
  )

  df_plot$site <- factor(df_plot$site,
                         levels = c("Ambrozio_2012", "Ambrozio_2019", "Caldart_2016", "This_work"),
                         labels = c("Ambrozio (2012)", "Ambrozio (2019)", "Caldart (2016)", "Este trabalho"))
g_ludico
