# Protokolle zum Grundpraktikum für Physikstudierende

**Beschreibung**

Aus [Moodle](https://moodle.tu-dortmund.de/):
> Zur Ausbildung eines Physik Studierenden an der [TU-Dortmund](https://www.tu-dortmund.de/) gehört im
> Bachelor Studiengang ein zweisemestriges Grundpraktikum, in dem an einfachen Standardversuchen unter
> anderem experimentelle Methoden in der Physik, Fehlerrechnung, das Schreiben von Protokollen sowie
> der Umgang mit physikalischen Geräten und Daten gelernt wird.

Die Struktur dieses Projekts und die grundlegende Methodik sind an den
[Toolbox-Workshop](https://toolbox.pep-dortmund.org/notes.html) von
[PeP et al. e.V.](https://pep-dortmund.org/) angelehnt. Als Hilfe stellt die
[Fachschaft](https://fachschaft-physik.tu-dortmund.de/wordpress/studium/praktikum/altprotokolle/)
Altprotokolle zur Verfügung.

**Autoren**

Fritz Agildere ([fritz.agildere@udo.edu](mailto:fritz.agildere@udo.edu)) und
Ben Brüggemann ([ben.brueggemann@udo.edu](mailto:ben.brueggemann@udo.edu))

**Struktur**

Die Protokolle werden mit `make` als PDF-Datei ausgegeben. Um alle Versuche auf einmal zu kompilieren,
kann `make` im Hauptverzeichnis ausgeführt werden. Hier wird auch die allgmeine Konfiguration vorgenommen.
Die Unterverzeichnisse übernehmen diese standardmäßig. Die einzelnen Versuche enthalten wiederum die
Verzeichnisse `build`, in dem sich alle generierten Dateien befinden, und `content`, das der Struktur
des Protokolls entspricht:

1. Theorie
2. Durchführung
3. Auswertung
4. Diskussion

Zur graphischen Darstellung und um abgeleitete Messwerte automatisch zu berechnen, werden `python` Skripte
mit den entsprechenden Bibliotheken genutzt. Die Dokumente werden unter Anwendung von `lualatex` gebaut.

Das Projekt *Anfängerpraktikum* ist mit GNU/Linux kompatibel.
