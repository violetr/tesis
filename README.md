### TESIS
código en R usado en la tesis de licenciatura :rocket:

## Código de los métodos utilizados

A continuacion damos los nombres de las funciones usadas. Los archivos se encuentran en la carpeta ``metodos``.


# Regresión lineal

* lasso, ridge: ``glmnet::glmnet/cv.glmnet``

En el archivo ``regresion_lineal.R</tt> se encuentran las siguientes funciones:

* adaptive lasso cv: ``cv.adaptive.lasso</tt>

* thresholded lasso cv: ``cv.tau.threslasso``


# Regresión logística


* regresion logistica con penalizacion: ``glmnet::glmnet/cv.glmnet ``



#Estimación del grafo normal


* GLasso: ``glasso::glasso``


En el archivo ``seleccion_modelo_normal.R`` se encuentran las siguientes funciones:


* Adaptive GLasso: ``adaptive.glasso``

* (Adaptive)GLasso CV: ``cv.glasso``

* Nodewise regression: ``nodewisereg``

* Stability Selection: ``stability.selection.grafico``

* Gelato (nodewise+threshold): ``cv.Gelato``


# Estimaci\'on del grafo discreto

* Chow-Liu: ``gRapHD::minForest``

* pseudo-likelihood: (http://github.com/yoshiomori/neighborhoods.git)


En el archivo ``seleccion_modelo_binario.R`` se encuentran las siguientes funciones:


* Nodewise Logistic Regression: ``nodewiselogreg``
(con criterio EBIC, CV, stability selection)


## Simulaciones y an\'alisis de datos

# Modelo normal

* ``3_2_1_3_tradeoff.R``: gr\'afico error, sesgo y varianza vs. $\lambda$ para el m\'etodo *lasso*.

* ``3_7_1_reg.R``: simulaci\'on regresi\'on lineal con *lasso*, *adaptive*, *thresholded* y *ridge*.


* ``3_7_2_consistencia.R``: gr\'afico de la consistencia en la estimaci\'on de la matriz de precisi\'on  para normal simulada.


* ``3_7_3_grafperdida.R``: gr\'afico de la p\'erdida con *GLasso``, *Adaptive GLasso* y *Gelato*``* para \cite{arabidopsis`` y normal simulada.

* ``3_7_4_ribo.R``: selecci\'on de modelo para ``hdi::riboflavin`` con *Stability Selection*. 



# Modelos discretos


* ``4_2_1_prostate.R``: clasificaci\'on de los datos ``spls::prostate``.

* ``4_4_1_palabras.R``: selecci\'on de modelo para datos \cite{news``.


* ``4_4_2_dipus.R``: selecci\'on de modelo para datos  de diputados en \cite{decadavotada``.


* ``4_4_3_bin.R``: simulaci\'on de selecci\'on de modelo para datos binarios.
