# Installation
git clone https://github.com/JoaquinMe/barque
# Dependencies
- GNU/Linux or OSX
- bash 4+
- python 3.5+
- python distutils package
- R 3+
- java
- gnu parallel
- flash (read merger) v1.2.11+
- vsearch v2.14.2+

# Preparación
- Instalar dependencias
- Descargar barque

# Preparar muestras 
./descargar_multiple.sh accessions.txt 04_data
python3 rename_script.py 04_data

# Primers
Editar 02_info/primers.csv para dar información de los primers
El default es el que vá. No hace falta que lo cambies

#BDs
La base de datos está bien

# Parámetros
Modificar 02_info/barque_config.sh de acuerdo a tu experimento

Correr barque:
./barque 02_info/barque_config.sh
