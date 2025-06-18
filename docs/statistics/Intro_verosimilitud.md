---
title: Introducción a la inferencia con la función de verosimilitud
layout: home
nav_order: 11
index: true
redirect: false
parent: Tutorials
math: katex
---

## Introducción


Ari Martínez pasa muchos veranos en Perú haciendo trabajo de campo en la jungla del Amazonas. En 2011 Ari observa un comportamiento interesante en las bandadas con especies mixtas. Si observaba a los pájaros que se alimentan en el suelo (llamados en inglés dead-leaf gleaning) los pájaros congelan su comportamiento por largo tiempo después de escuchar una señal de alarma. Mientras tanto los pájaros mosqueros (flycatchers) que se alimentan en las copas de los árboles ni siquiera se preocupaban de la señal de alarma. 

Ari empezó a medir el tiempo de este comportamiento. Midiendo en segundos cuánto tiempo los pájaros se congelaban. En su experimento, él sonaba la alarma y empezaba el conteo de los segundos, el conteo acaba cuando los pájaros seguían comiendo. 

La muestra de Ari es lo que vamos a usar para entender conceptos básicos de probabilidades y la función de verosimilitud. 

## Descarga los datos

 + Lista de tiempos y tipo de pájaros [aquí](files/birdalarms.csv)

``` 
#Cuál es tu directorio de trabajo? O puedes comenzar con un nuevo proyecto .Rproj
#setwd()
bird.alarms<-read.csv("birdalarms.csv",stringsAsFactors=FALSE)
head(bird.alarms)
```

### Variables aleatorias


### El modelo

A la variable aleatoria del tiempo de espera le vamos a asignar una función de densidad de probabilidad llamada **distribución exponencial**. 

La distribución exponencial tiene un parámetro que es una tasa llamada $$1/\theta$$. La $$\theta$$ se interpreta como el promedio (o la esperanza, en terminos formales estadísticos) del tiempo que tardan los pájaros en descongelarse. 


### La función de la verosimilitud

La función de verosimilitud, o verosimilitud en términos sencillos es
$$L(\theta; \mathbb{X}) = P(\mathbb{X}|Modelo (\theta))$$

y se interpreta como **la probabilidad de la muestra dado el modelo o la hipótesis**. La función de verosimilitud tiene dos partes
1. El modelo: Representado por la distribución exponencial con parámetro $$1/\theta$$.
2. Los datos: $$\mathbb{X}$ esto representa la muestra de nuestras observaciones. En este ejemplo la muestra es el tiempo que tardo cada uno de nuestros pájaros en descongelarse.

La verosimilitud es una función sencilla de entender, lo que nos dice es que estamos calculando la probabilidad de observar todos estos comportamientos de los pájaros bajo un modelo razonable. Nosotros podemos cambiar el modelo cómo queramos, pero los datos es lo que hemos observado.

Calculemos esta verosimilitud en R

``` 
likelihood.function<- function(parameter, observations){
probabilities <-dexp(observations, rate=1/parameter,log=FALSE)
L <- prod(probabilities)
return(-L)
}
```

En esta función
+ Qué ponemos como parámetros?
+ Qué resulta y por qué lo que resulta es negativo?


### Un mini ejemplo con dos observaciones

Tenemos una muestra con $$n=$$, la muestra es $\mathbb{X}={2, 10}$$ . Antes de avanzar pensemos, cuál sería el número promedio de segundos que tengo que esperar para que los pájaros de esta muestra se descongelen?

``` 
mle<-?
```

Con una muestr más grande o modelos más complicados la verosimilitud se complica.
Veamos como optimizar

``` 
(likelihood.optimization<-optimize(f=likelihood.function, interval=c(0,20), observations=c(2,10)))
```


La función `optimize` se utiliza cuando queremos optimizar solo un parámetro. Si mi modelo tiene más de un parámetro se utiliza `optim` o
`nloptr` (del paquete nloptr)

Veamos lo que resulta

+ Valor máximo verosímil (MLE):
+ Valor de la función de verosimilitud evaluado en el MLE:

### El ejemplo completo de Ari!


Para los pájaros que se alimentan en el suelo (dead-leaf gleaning) tenemos que sus datos son

``` 
dl.forager<- bird.alarms$Response[which(bird.alarms$Forage=="DL")]
```

y para los mosqueros (flycatchers) tenemos que 

```
f.forager<-bird.alarms$Response[which(bird.alarms$Forage=="F")] #11
```

#### Calcula la optimización para estos ejemplos

Dead-leaf gleaning:

+ Valor máximo verosímil (MLE):
+ Valor de la función de verosimilitud evaluado en el MLE:

Flycatchers:
+ Valor máximo verosímil (MLE):
+ Valor de la función de verosimilitud evaluado en el MLE:

Qué observas aquí? Se pueden comparar estos valores?


### Evidencia de que estos grupos de pájaros se comportan distinto

Las funciones de verosimilitud no solo son su valor en el máximo verosímil, sino que representan plausibilidad. Esto significa que nos permiten calcular la plausibilidad de diferentes valores del parámetro bajo nuestro modelo y entender si los datos tienen buena probabilidad bajo diferentes valores del parámetro. Por eso es importante explorar la verosimilitud más.

Exploremos entre (0, 50) segundos de congelamiento

``` 
parameter.vals<-seq(0.0001,50,0.01) #Intervalo
long<-length(parameter.vals) # largo del intervalo
# Evaluando la verosimilitud para cada uno de esos segundos
p.likelihoodf<-rep(0,long)
for (i in 1:long){
p.likelihoodf[i]<- -likelihood.function(parameter.vals[i],observations=f.forager) # Remeber is negative so we need to add a sign
}
```

Verosimilitud en diferentes escalas

``` 
par(mfrow=c(1,3))
# Directo como la calculamos

plot(parameter.vals, p.likelihoodf, type="l",main="Likelihood for flycatchers",xlab=expression(theta),ylab="Likelihood", lwd=2,xlim=c(0,5))


# En escala logaritmica
plot(parameter.vals, log(p.likelihoodf), type="l",main="log-likelihood for flycatchers",xlab=expression(theta),ylab="Likelihood",lwd=2,xlim=c(0,5))

# Verosimilitud relative: Verosimilitud dividida por el valor maximo en el MLE
plot(parameter.vals, p.likelihoodf/likelihoodval.f, type="l",main="Relative likelihood for flycatchers",xlab=expression(theta),ylab="Likelihood",lwd=2,xlim=c(0,5))
```

![](/likelihoodimages/likelihoodforf-1.png)


#### Ahora para los que se alimentan en el suelo
``` r
p.likelihooddl<-rep(0,long)
for (i in 1:long){
p.likelihooddl[i]<- -likelihood.function(parameter.vals[i],observations=dl.forager)
}

par(mfrow=c(1,3))
plot(parameter.vals, p.likelihooddl, type="l",main="Likelihood for flycatchers",xlab="Rate Parameter",ylab="Likelihood",lty=2,col="red",lwd=2)

plot(parameter.vals, log(p.likelihooddl), type="l",main="log-likelihood for flycatchers",xlab="Rate Parameter",ylab="Likelihood",lty=2,col="red",lwd=2)

plot(parameter.vals, p.likelihooddl/likelihoodval.dl, type="l",main="Relative likelihood for flycatchers",xlab="Rate Parameter",ylab="Likelihood",lty=2,col="red",lwd=2)
```


Se ve como esto

![](/likelihoodimages/unnamed-chunk-3-1.png)

### Los grupos son distintos?- Test the razones de verosimilitud (Likelihood ratio test)


Usualmente, en estadística frecuentista alguien hubiera sugerido utilizar un "t-test" y probar que estos grupos son diferentes. Pero en un t-test la muestra se asume que viene de otra distribución distinta (la distribución Normal) y nosotros asumimos otra distribución y tenemos muy poquitos datos. Entonces, cómo pruebo que la distribución es distinta? La verosimilitud tiene muy buenas propiedades, utilizemosla. 

``` 
plot(parameter.vals, p.likelihoodf/likelihoodval.f, type="l",main="Evidence for responses",xlab="Rate Parameter",ylab="Likelihood")
lines(parameter.vals, p.likelihooddl/likelihoodval.dl,lty=2,col="red",lwd=2)
legend(x=35,y=0.8, col=c("black","red"),legend=c("flycatcher","dead-leaf"),lty=1:2)
```

![](/likelihoodimages/evidenceplots-1.png)

### Intervalos de verosimilitud-confianza

De la teoría estadistica: Un resultado asintótico de la verosimilitud relativa. Es que podemos cortar la verosimilitud relativa en el 14.6% y obtenemos intervalos de verosimilitud-confianza del 95% y son más acertados que cualquier otro intervalo. 

``` 
plot(parameter.vals, p.likelihoodf/likelihoodval.f, type="l",main="Evidence for responses",xlab="Rate Parameter",ylab="Likelihood")
lines(parameter.vals, p.likelihooddl/likelihoodval.dl,lty=2,col="red",lwd=2)
legend(x=35,y=0.8, col=c("black","red"),legend=c("flycatcher","dead-leaf"),lty=1:2)
abline(h=0.146,lty=3)
```

![](/likelihoodimages/likelihoodintervals-1.png)

### Test de razones de verosimilitud (likelihood-ratio test)
Si todos nuestro revisores entendieran verosimilitud, estarían satisfechos de ver nuestros resultados. Sin embargo, muchos esperan el famoso p-valor. Vamos a hacerlo entonces utilizando las propiedades de verosimilitud. Utilizaremos un estadistico de prueba llamado **LRT**.

$$LRT = -2(log (L (H_0))- log(L(H_1)))

Para hipotesis nula $$H_0$$  vamos a asumir que todos los pajaros tienen el mismo comportamiento. 

```
mle.H0<-mean(c(dl.forager,f.forager))
loglike.H0<- log(-likelihood.function(parameter=mle.H0,observations=c(dl.forager,f.forager)))
```

Y para la alternativa $$H_1$$ asumimos que cada grupo tiene distinto comportamiento.
``` r
(loglike.H1<-log(likelihoodval.f*likelihoodval.dl))
 -77.5975
```

Estos dos modelos son anidados!

Ahora calculamos el estadístico de prueba
``` 
(LRT<- -2*(loglike.H0-loglike.H1))
```

Que tiene una distribución chi-square con un solo grado de libertad. El p-valor de la prueba es

``` 
(p.value<-pchisq(LRT,df=1,lower.tail=FALSE))
```



1. Martínez, A.E. and Zenil, R.T., 2012. Foraging guild influences dependence on heterospecific alarm calls in Amazonian bird flocks. Behavioral Ecology, 23(3), pp.544-550.

2. Strug, L.J., 2018. The evidential statistical paradigm in genetics. Genetic epidemiology, 42(7), pp.590-607.


