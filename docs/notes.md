# Notebook1

Carica i dati e crea le directory

Seleziona un "target spectra": usa il metodo SAMOS.read_SAMI_mosaic(file) per creare un hdu, ed i suoi dati li salva come "Target_sectra", poi usa il metodo SAMOS.CR_correct(Spectra) su di essi

Crea l'oggetto "flat", corregge per CR e lo mostra (luce grossa a destra)

Stessa procedura per i "bias" e poi somma i diversi file di bias

Poi fa flat - bias (a cosa serve?) e dopo spectra - bias (sembra non andare bene). Fa una procedura in più per correggere spectra con il suo proprio bias (che cos'è quel plot?)

Ora carica ARC LAMP e la plotta, idem per una ARC con il suo Dark --> l'immagine è buona (si vede bene all'occhio) --> allora somma su 1400 bin (in orizzontale) e viene il plot con le posizioni delle slitte

Poi fa un diff per trovare la posizione giusta... c'è il problema delle slits troppo vicine, fa due scamotagge per definire la loro posizione

Poi fa il plot delle slits

Dopo dei flats --> perchè?

Dopo del flat bias blur --> perchè?

Dopo un loop su tutte le slits --> MASKS



### Commenti:

target_mode e File_type --> devono essere letteralmente così? elenco dei possibili valori all'inizio!


### Domande:

cos'è quella luce grossa a destra nel flat?

Cosa rappresenta quel plot del flat? perchè quel picco a 30000?



# Notebook2

### Inputs

### Outputs

### Function
