### TESIS
código en R usado en la tesis de licenciatura :rocket:

## C\`odigo de los m\`etodos utilizados

A continucaci'on damos los nombres de las funciones usadas. Los archivos se encuentran en la carpeta <\tt>metodos<\tt>.


# Regresi\'on lineal

* lasso, ridge: <tt>glmnet::glmnet/cv.glmnet<\tt>

En el archivo <tt>regresion_lineal.R</tt> se encuentran las siguientes funciones:

* adaptive lasso cv: <tt>cv.adaptive.lasso</tt>

* thresholded lasso cv: <tt>cv.tau.threslasso<\tt>


# Regresi\'on log\'istica


* regresion logistica con penalizacion: <tt>glmnet::glmnet/cv.glmnet <\tt>



#Estimaci\'on del grafo normal


* GLasso: <tt>glasso::glasso<\tt>


En el archivo <tt>seleccion_modelo_normal.R<\tt> se encuentran las siguientes funciones:


* Adaptive GLasso: <tt>adaptive.glasso<\tt>

* (Adaptive)GLasso CV: <tt>cv.glasso<\tt>

* Nodewise regression: <tt>nodewisereg<\tt>

* Stability Selection: <tt>stability.selection.grafico<\tt>

* Gelato (nodewise+threshold): <tt>cv.Gelato<\tt>


# Estimaci\'on del grafo discreto

* Chow-Liu: <tt>gRapHD::minForest<\tt>

* pseudo-likelihood: (http://github.com/yoshiomori/neighborhoods.git)


En el archivo <tt>seleccion_modelo_binario.R<\tt> se encuentran las siguientes funciones:


* Nodewise Logistic Regression: <tt>nodewiselogreg<\tt>
(con criterio EBIC, CV, stability selection)


## Simulaciones y an\'alisis de datos

# Modelo normal

* <tt>3_2_1_3_tradeoff.R<\tt>: gr\'afico error, sesgo y varianza vs. $\lambda$ para el m\'etodo *lasso*.

* <tt>3_7_1_reg.R<\tt>: simulaci\'on regresi\'on lineal con *lasso*, *adaptive*, *thresholded* y *ridge*.


* <tt>3_7_2_consistencia.R<\tt>: gr\'afico de la consistencia en la estimaci\'on de la matriz de precisi\'on  para normal simulada.


* <tt>3_7_3_grafperdida.R<\tt>: gr\'afico de la p\'erdida con *GLasso<\tt>, *Adaptive GLasso* y *Gelato*<\tt>* para \cite{arabidopsis<\tt> y normal simulada.

* <tt>3_7_4_ribo.R<\tt>: selecci\'on de modelo para <tt>hdi::riboflavin<\tt> con *Stability Selection*. 



# Modelos discretos


* <tt>4_2_1_prostate.R<\tt>: clasificaci\'on de los datos <tt>spls::prostate<\tt>.

* <tt>4_4_1_palabras.R<\tt>: selecci\'on de modelo para datos \cite{news<\tt>.


* <tt>4_4_2_dipus.R<\tt>: selecci\'on de modelo para datos  de diputados en \cite{decadavotada<\tt>.


* <tt>4_4_3_bin.R<\tt>: simulaci\'on de selecci\'on de modelo para datos binarios.
