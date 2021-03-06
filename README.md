# TESIS
código en R usado en la tesis de licenciatura :rocket:

## Código de los métodos utilizados

A continuacion damos los nombres de tolas las funciones usadas. Los archivos de código se encuentran en la carpeta ``metodos``.


### Regresión lineal

* lasso, ridge: ``glmnet::glmnet/cv.glmnet``

* stability selection: ``hdi::stability``

En el archivo ``regresion_lineal.R`` se encuentran las siguientes funciones:

* adaptive lasso cv: ``cv.adaptive.lasso``

* thresholded lasso cv: ``cv.tau.threslasso``


### Regresión logística


* regresion logistica con penalizacion: ``glmnet::glmnet/cv.glmnet ``



### Estimación del grafo normal


* GLasso: ``glasso::glasso``


En el archivo ``seleccion_modelo_normal.R`` se encuentran las siguientes funciones:


* Adaptive GLasso: ``adaptive.glasso``

* (Adaptive)GLasso CV: ``cv.glasso``

* Nodewise regression: ``nodewisereg``

* Stability Selection: ``stability.selection.grafico``

* Gelato (nodewise+threshold): ``cv.Gelato``


### Estimación del grafo discreto

* Chow-Liu: ``gRapHD::minForest``

* pseudo-likelihood: (http://github.com/yoshiomori/neighborhoods.git)


En el archivo ``seleccion_modelo_binario.R`` se encuentran las siguientes funciones:


* Nodewise Logistic Regression: ``nodewiselogreg``
(con criterio EBIC, CV, stability selection)


## Simulaciones y análisis de datos

Los siguientes archivos se encuentran en la carpeta ``analisis``.

### Modelo normal

* ``3_2_1_3_tradeoff.R``: gráfico error, sesgo y varianza vs. lambda para el método *lasso*.

* ``3_7_1_reg.R``: simulación regresión lineal con *lasso*, *adaptive*, *thresholded* y *ridge*.


* ``3_7_2_consistencia.R``: gráfico de la consistencia en la estimación de la matriz de precisión  para normal simulada.


* ``3_7_3_grafperdida.R``: gráfico de la pérdida con *GLasso*, *Adaptive GLasso* y *Gelato* para ``arabidopsis`` y normal simulada.

* ``3_7_4_ribo.R``: selección de modelo para ``hdi::riboflavin`` con *Stability Selection*. 



### Modelos discretos


* ``4_2_1_prostate.R``: clasificación de los datos ``spls::prostate``.

* ``4_4_1_palabras.R``: selección de modelo para datos ``news``.

* ``4_4_2_dipus.R``: selección de modelo para datos  de diputados en ``decadavotada``.

* ``4_4_3_bin.R``: simulación de selección de modelo para datos binarios.
