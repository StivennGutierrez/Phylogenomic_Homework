## ANÁLISIS FILOGENÓMICO DE PROTEOMA DE 23 SPP DE VERTEBRADOS

### Resumen del workflow de obtención de árboles de coalescencia y de máxima verosimilitud a partir de proteomas de especies de vertebrados:

En un inicio, claro esta lo principal es obtener los datos potencialmente analizables por filogenómica (en este caso correspondientes a proteomas de 23 spp de vertebrados:elephant:). Luego, al ya tener nuestras secuencias en una carpeta de trabajo, activamos el ambiente `conda` :gear: :desert:. Lo primero es iniciar con la inferencia de ortólogos, nuestro objetivo es rastrear regiones ortólogas entre proteínas. Para ello empleamos el software __Orthofinder__ con el comando `orthofinder` especificandole la carpeta :card_index_dividers: con los respectivos proteomas :dna:, nos generará de esta forma un archivo de secuencia para cada uno en formato __.tsv__. Por último y para facilitar la concatenación, empleamos expresiones regulares para cambiar la estructura de los nombres de las secuencias a `Género_especie` lo más conveniente es usar `sed`. Después, es importante hacer un filtraje de calidad🔎 con `PREQUAL` para descartar todos los posibles errores de ensamblaje o anotación. Con ello, para corroborar homología entre aa realizamos un __alineamiento múltiple__ por `mafft`. Además, removemos:scissors: posibles gaps generados con `BMGE`:wastebasket:. Por útlimo, realizamos un alinamiento cancatenando:abacus: por medio de __FASconCAT__ y `perl`, para así obtener 23 taxa con 21 genes. Con esa información generamos dos árboles por máxima verosimilitud: 1. un solo locus que evoluciona de la misma forma:round_pushpin: y 2. por partición según ubicación de genes, usando `IQTREE`:bar_chart:. Así como un árbol de coalescencia usando `ASTRAL`:comet:.

#### 1) Árbol de máxima verosimilitud (no segmentado):
![unpartitioned treefile](https://github.com/StivennGutierrez/parcial_bioinformatica/assets/128840301/65473ddb-c357-4b59-a851-ab12b233e02a)

----

### Resultado:
En general, el árbol señala la formación de distintos clados como el que incluye las especies *Oreochromis niloticus*:tropical_fish:, *Takifugu rubripes*:fish:, *Danio rerio*:fried_shrimp: y *Lepisosteus oculatus*:shrimp:, donde la especie *T. rubripes* señala mayor acumulación de cambios evolutivos a tráves del tiempo :hourglass_flowing_sand:. Cabe resaltar que este grupo monofilético muestra apariciones y cambios evolutivos más lejanos entre cada especie.

Asimismo, es posible identificar otros clados con apariciones y cambios evolutivos no muy lejanos entre sí: el primero, conformado por las especies *Pelodiscus sinensis*:turtle:, *Anolis carolinensis*:lizard:, *Taeniopygia guttata*:bird:, *Meleagris gallopavo*:turkey: y *Gallus gallus*	:rooster:; el segundo, conformado por *Homo sapiens*:tipping_hand_person:, *Mus musculus*:mouse2:, *Canis lupus familiaris*:dog2:, *Dasypus novemcintus*:hedgehog: y *Loxodonta africana*:elephant:, siendo este último grupo monofilético más reciente.

Pero si nos fijamos bien, podemos observar que las especies *Neoceratodus forsteri*:penguin:, *Protopterus annectens*:shrimp: y *Lepidosiren paradoxa*:shrimp: se separan de los clados anteriores un tiempo evolutivo considerable, siendo *N. forsteri* la especie más reciente :new:. Además, la especie *Callorhinchus milli*:shark: sería la especie que menos características compartiría con el resto de clados, sin embargo debe presentar una característica que lo separó del primer clado mencionado, y que a lo largo del tiempo ha permanecido en un constante cambio evolutivo:electron:.

---

#### 2) Árbol de máxima verosimilitud (segmentado):
![partitioned treefile](https://github.com/StivennGutierrez/parcial_bioinformatica/assets/128840301/96ca351e-4e38-476b-972a-caff3a851152)

----

### Resultado:
En este caso, los valores de bootstrap nos permiten corroborar la fiabilidad de las relaciones filogenéticas obtenidas. De esta forma, atendiendo a lo mencionado en el análisis anterior, la ubicación de *A. carolinensis* realmente es incierta teniendo en cuenta el valor de soporte tan bajo (44.6). De resto, todas las relaciones filogenéticas señalan valores de bootstrap bastante confiables. Y haciendo un barrido de las distancias de rama, podemos observar que al generar particiones dentro del análisis, la sensibilidad frente a los cambios evolutivos es mayor, ya que posiblemente se pudieron contemplar cambios en las secuencias dentro de un modelo de evolución para distintos locus.

---

#### 3) Árbol de coalescencia:
![species_tree_ASTRAL tre](https://github.com/StivennGutierrez/parcial_bioinformatica/assets/128840301/fe5062a5-2389-4ae2-8b73-322aa6dac1cb)

----

### Resultado:
Ya para un análisis de coalescencia, el árbol nos indica los tiempos de coalescencia y las longitudes de rama para el caso de cada especie según la señal de cambio evolutivo en unidades coalescentes. De este modo, es posible observar que las inferencias de las relaciones filogenéticas son diferentes a las propuestas por el algoritmo de máxima verosimilitud. De hech, para la primera división de clados se evidencia un claro evento de politomía, que teniendo en cuenta el valor bajo de bootstrap, solo no se pudo inferir una clara filiación de ese ancestro frente al resto, es decir, aún no es claro el tiempo de coalescencia en el que los ancestros de los respectivos clados dieron inició a estos clados.

---

#### Diferencias entre los tres tipos de árboles:
| **Aspecto** | **Árbol de máxima verosimilitud sin particiones** | **Árbol de máxima verosimilitud particionado** | **Árbol de coalescencia** |
| --- | --- | --- | --- |
| Valor de soporte | Bootstrap | Bootstrap | Basado en probabilidades posteriores o valores de bootstrap por medio de múltiples árboles generados por estocasticidad de los procesos coalescentes |
| Representación rama | Se basa en las distancias filogenéticas | Se basa en las distancias filogenéticas | Representa unidades de tiempo coalescente en generaciones o años |
| Interpretación y resultados | Representa el árbol que mejor ubica las distancias evolutivas entre especies | Representa la estimación más adecuada para representar la historia evolutiva dadas las distancias evolutivas de cada especie | Representa la historia genealógica más probable según los datos de secuencia y las inferencias propuestas por un modelo coalescente |
| Uso | Estimar relaciones filogenéticas para secuencias lejanamente relacionadas | Estimar relaciones filogenéticas para secuencias estrechamente relacionadas | Establecer relaciones filogenéticas a nivel poblacional a partir de eventos de coalescencia |
