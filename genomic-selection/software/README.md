## Software de flujo de trabajo para Selección Genómica
El software realiza el análisis comparativo para 14 algoritmos de selección genómica. Trabaja en línea de comandos y recibe como parámetros de entrada tres archivos:
- Primero, Los genotipos de la población de entrenamiento
- Segundo, Los genotipos de la población de prueba
- Tercero, Los fenotipos de la población de entrenamiento

El software produce por cada análisis dos tipos de archivos, uno texto en formato CSV con la tabla con de las comparaciones y un PDF con la gráfica comparativa.

El llamado del software desde la línea de comandos es el siguiento
```
 GSPipeline.R <Training genotype> <Testing genotype> <Training phenotype>
```
