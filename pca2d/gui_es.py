"""The window in Spanish (added 2026-09-16), beside gui.EN and gui.FR.

Every key of gui.EN, with the same placeholders in the same order. Written in
the formal register (usted), in the Spanish spoken where the flag on the
language button is from.
"""

ES = {
    "subtitle": "preclean en dos marcos, luego LBL. Elija una raíz de datos,"
                " elija objetos, revise el comando y ejecútelo.",
    "data_dir": "raíz de datos",
    "config": "configuración",
    "fits_dir": "disco de productos (opcional)",
    "help_fits_dir":
        "Dónde se GUARDAN los productos de una ejecución: cada carpeta de"
        " ejecución pasa a ser un enlace a este disco, de modo que los"
        " espectros corregidos, el ajuste y las figuras quedan allí y no en el"
        " disco interno. Vacío, todo queda bajo la raíz de salida, lo cual"
        " funciona. Una ruta existe en UNA sola máquina, así que este campo se"
        " vacía cada vez que lo que nombra no está, y una ejecución cuyo disco"
        " falta se detiene en lugar de llenar sin avisar el disco interno.",
    "log_no_disk":
        "la configuración nombra %s como el disco donde guardar los productos,"
        " y no está: el campo queda vacío, así que esta ejecución guardaría"
        " todo bajo la raíz de salida. Rellénelo si el disco debería estar"
        " montado.",
    "out_dir": "raíz de salida (opcional)",
    "lbl_dir": "carpeta de salida de LBL",
    "help_lbl_dir":
        "El árbol propio de LBL (su DATA_DIR): las carpetas science que lee, y"
        " las plantillas, máscaras, tablas por línea y rdb que escribe. Se"
        " propone como lbl bajo la raíz de salida, y sigue a esa raíz hasta que"
        " usted escriba otra. Un solo árbol para todas las ejecuciones bajo una"
        " raíz, ya que el LBL del objeto entregado es el mismo para todas y"
        " toma horas; apúntelo a un árbol existente y ese árbol se usa tal"
        " cual.",
    "browse": "Examinar",
    "rescan": "Releer",
    "objects": "objetos",
    "settings": "ajustes",
    "stages": "etapas",
    "command": "el comando que se ejecutará",
    "output": "salida",
    "col_object": "objeto",
    "col_files": "archivos",
    "col_instrument": "instrumento",
    "col_snr": "S/R",
    "col_exptime": "exp. (s)",
    "col_mag": "mag",
    "run_name": "nombre de la reducción",
    "auto": "auto",
    "help_auto_button":
        "Vuelve a proponer un nombre a partir de los objetivos y de los ajustes"
        " tal como están: los objetivos, y luego seis caracteres de un hash de"
        " todo lo que hace de esta reducción un resultado distinto. Los mismos"
        " parámetros dan los mismos seis; un número distinto da otros seis.",
    "timeline": "cuándo se observaron las campañas marcadas",
    "timeline_none": "marque un objetivo para ver cuándo se observó",
    "timeline_note": "se conserva de %s a %s: %d exposiciones de %d",
    "whole": "todo",
    "min_rjd": "desde",
    "max_rjd": "hasta",
    "exists": "esta ejecución ya existe",
    "help_run_name":
        "Un nombre para esta ejecución. Sus productos van a"
        " <raíz de salida>/_NOMBRE/ y su objeto LBL es"
        " <objeto>_PCA2D_<M-N>_NOMBRE, de modo que dos ejecuciones de los"
        " mismos objetivos con ajustes distintos nunca escriben en una misma"
        " carpeta ni bajo un mismo nombre LBL, donde LBL mediría la mezcla sin"
        " decir nada. Vacío es la ruta nominal. Se propone a partir de las"
        " fechas cuando están fijadas, y puede ser cualquier cosa.",
    "help_dates":
        "Conserva solo las exposiciones entre estas dos fechas, en fecha"
        " juliana reducida (BJD - 2400000). Lo excluido no se ajusta ni se"
        " corrige. Así se recorta una campaña a una temporada, cosa que a veces"
        " pide la cobertura baricéntrica: dos cortes de 14 noches de TOI-4552"
        " con la misma señal a ruido difieren en un factor sesenta en"
        " cobertura, y la corrección cambia de signo entre ellos.",
    "berv": "cobertura baricéntrica de los objetivos marcados",
    "berv_none": "marque un objetivo para ver su cobertura baricéntrica",
    "berv_wait": "aún no se ha leído ninguna velocidad baricéntrica",
    "berv_building": "EN CONSTRUCCIÓN: todavía leyendo",
    "berv_note":
        "cobertura efectiva %.0f km/s, extensión %.0f, como máximo %s posible"
        " para la latitud eclíptica de este objetivo; intervalos de %.0f km/s,"
        " %d exposiciones",
    "help_berv":
        "Qué parte del intervalo baricéntrico cubren realmente las campañas"
        " marcadas, en intervalos de 3 km/s, un color por instrumento y"
        " apilados. Es la cantidad que decide si la corrección ayuda: el bloque"
        " del observador solo es identificable porque la estrella se desplaza"
        " a través de él. Dos cortes de 14 noches de TOI-4552 con la misma"
        " señal a ruido, uno de 42.8 km/s de extensión y el otro de 0.7,"
        " pasaron de una ganancia de un factor dos a una pérdida de un factor"
        " 2.4. La cobertura EFECTIVA cuenta solo los intervalos que contienen"
        " una exposición, de modo que a una campaña observada en dos extremos y"
        " en ningún punto intermedio no se le acredita el intervalo entre"
        " ellos. El tercer número es lo que el cielo PERMITE: el |BERV| de un"
        " objetivo nunca supera 29.78 cos(latitud eclíptica) km/s, así que"
        " TOI-1452, a +80.5 grados, solo puede abarcar 9.8 km/s y ninguna"
        " cantidad de observaciones separará los dos marcos para él. Una"
        " cobertura menor que la extensión indica huecos, que más noches pueden"
        " llenar; una extensión menor que lo posible indica una campaña joven;"
        " un posible pequeño indica que el objetivo no es el adecuado.",
    "run": "Ejecutar",
    "stop": "Detener",
    "dry": "Simulación",
    "export": "Exportar YAML...",
    "savelog": "Guardar registro...",
    "openout": "Abrir salidas",
    "openpdf": "Abrir PDF de compilación",
    "savedefaults": "Guardar como predeterminados...",
    "all": "todos",
    "none": "ninguno",
    "idle": "en espera",
    "running": "en curso",
    "lang": "Español",
    "snr_berv": "señal a ruido frente a la velocidad baricéntrica",
    "snr_none": "marque un objetivo para ver dónde están sus mejores noches",
    "help_snr_berv":
        "Un punto por espectro, en el color de su estrella: la señal a ruido"
        " que midió APERO frente a la velocidad baricéntrica en que se tomó. El"
        " histograma de arriba dice qué velocidades cubre una campaña; este"
        " dice CON QUÉ las cubre. El ajuste pondera un espectro por"
        " 1/sigma^2, así que un intervalo cubierto en un extremo por las peores"
        " noches de una campaña no es el intervalo que el ajuste ve de verdad,"
        " y los dos marcos se separan peor de lo que la cobertura promete.",
    "tab_targets": "  objetivos  ",
    "tab_settings": "  ajustes  ",
    "tab_lbl": "  LBL  ",
    "tab_run": "  análisis  ",
    "tab_clean": "  limpieza  ",
    "tab_runs": "  pasadas  ",
    "runs_title": "todas las pasadas que contiene esta carpeta de salida",
    "help_runs":
        "Lo que se ha probado, leído de lo que cada pasada dejó escrito: una"
        " pasada registra los ajustes que resolvió, la orden que recibió y la"
        " hora en que empezó, así que esta lista son las pasadas mismas y no"
        " un historial que guarde la ventana. Elija una para ver todas sus"
        " opciones y abrir su PDF de compilación.",
    "runs_rescan": "Volver a leer",
    "help_runs_rescan":
        "Recorrer de nuevo la carpeta de salida y listar todas las pasadas que"
        " hay en ella. Una pasada en curso aparece en cuanto ha escrito su"
        " configuración resuelta, es decir, antes de su primera etapa.",
    "runs_open_pdf": "Abrir su PDF de compilación",
    "help_runs_open_pdf":
        "Abrir el único PDF de la pasada elegida en la lista: sus figuras, sus"
        " números y, una vez medido por LBL, sus velocidades. Lo escribe la"
        " etapa de figuras y de nuevo después de LBL.",
    "runs_open_folder": "Abrir su carpeta",
    "help_runs_open_folder":
        "Mostrar la carpeta de la pasada elegida en el explorador de archivos:"
        " la configuración resuelta, el ajuste, los espectros corregidos y el"
        " informe están todos ahí.",
    "runs_scanning": "leyendo la carpeta de salida...",
    "runs_found": "%d pasada(s) en %s",
    "runs_none": "aún no hay ninguna pasada en %s: una pasada escribe ahí su"
                 " configuración resuelta antes de su primera etapa",
    "runs_pick": "elija primero una pasada en la lista",
    "runs_no_pdf": "aún no hay PDF de compilación en %s",
    "runs_command": "la orden que recibió",
    "col_hash": "hash", "col_started": "empezó",
    "col_targets": "objetivos", "col_tag": "contadores",
    "col_components": "estrella + observador", "col_report": "PDF",
    "col_setting": "ajuste", "col_value": "esta pasada",
    "clean_title": "lo que pca2d ha dejado en estos discos",
    "help_clean":
        "Cada lugar donde este programa escribe bytes, medido. La caché de"
        " cubos, los temporales de un ajuste demasiado grande para la memoria"
        " y los registros de LBL pueden borrarse en cualquier momento: se"
        " rehacen a partir de lo que sigue aquí, y un botón los vacía todos."
        " Todo lo demás también puede borrarse, fila por fila: elíjala en la"
        " lista y Borrar la selección. Los ajustes, las plantillas, las"
        " máscaras y las tablas por línea solo vuelven ejecutando las horas que"
        " los produjeron, y los espectros corregidos, los informes y las"
        " velocidades solo volviendo a ejecutar todo desde los espectros; la"
        " ventana dice cuáles de esas horas está a punto de gastar antes de"
        " borrar nada. Una carpeta trasladada a otro disco se cuenta donde está"
        " realmente; un enlace a los espectros de otra persona no pesa nada,"
        " porque borrarlo no libera nada.",
    "clean_measure": "Medir",
    "clean_purge": "Borrar lo que puede irse",
    "clean_selected": "Borrar la selección",
    "help_clean_measure":
        "Recorre cada carpeta y suma lo que contiene. En un disco de red esto"
        " tarda unos segundos y se hace aparte, de modo que la ventana sigue"
        " utilizable mientras cuenta.",
    "help_clean_purge":
        "Vacía las carpetas marcadas como temporales o reconstruibles, y nada"
        " más: las que solo cuestan tiempo recuperar. Pregunta antes, e indica"
        " el total que va a liberar. Las carpetas mismas se quedan: la próxima"
        " ejecución espera encontrarlas. Para todo lo demás, elija las filas y"
        " use Borrar la selección.",
    "help_clean_selected":
        "Borra las filas elegidas en la lista, sean del tipo que sean: una"
        " caché de cubos, las plantillas que tomaron una tarde, las velocidades"
        " mismas. Mayúsculas o comando permite elegir varias. Nombra cada"
        " carpeta y dice lo que costaría recuperarla antes de borrar nada, y"
        " eso es lo único que hay entre usted y un disco con espacio libre.",
    "clean_measuring": "midiendo…",
    "clean_totals": "%s en total, de los cuales se pueden liberar %s",
    "clean_none": "nada medido todavía: pulse Medir",
    "clean_confirm_title": "borrar los archivos reconstruibles",
    "clean_confirm":
        "A punto de liberar %s en %d lugares:\n\n%s\n\nNada de esto es un"
        " resultado: los cubos se vuelven a leer desde los espectros, y LBL"
        " vuelve a escribir sus registros. ¿Continuar?",
    "clean_go": "Barrerlo",
    "clean_keep": "Conservarlo",
    "clean_freed": "%s liberados en %d lugares",
    "clean_nothing": "nada que liberar: aquí no hay temporales ni caché",
    "clean_pick":
        "elija primero una fila o varias en la lista: el botón borra lo"
        " elegido",
    "clean_while_running":
        "HAY UNA EJECUCIÓN EN CURSO, y lee de estas carpetas. Reconstruye lo"
        " que encuentre ausente, así que borrar ahora le cuesta a esa ejecución"
        " el tiempo de rehacerlo, en medio de la etapa en la que está.",
    "clean_confirm_pick_title": "borrar la selección",
    "clean_confirm_pick":
        "A punto de borrar %s en %d lugares:\n\n%s\n\n%s\n\n¿Continuar?",
    "clean_cost_rebuildable":
        "Todo se rehace a partir de lo que queda aquí: minutos, y una ejecución"
        " que vuelve a leer los espectros.",
    "clean_cost_expensive":
        "Una parte solo vuelve ejecutando lo que la produjo: un ajuste, o una"
        " pasada de LBL sobre cada exposición. Horas, no minutos.",
    "clean_cost_results":
        "UNA PARTE ES UN RESULTADO: espectros corregidos, un informe, o las"
        " velocidades que guarda un rdb. Nada aquí los rehace. Solo volver a"
        " ejecutar todo, desde los espectros, y eso son las horas que tomó la"
        " primera vez.",
    "col_size": "tamaño",
    "col_nfiles": "archivos",
    "col_kind": "tipo",
    "kind_scratch": "temporal",
    "kind_rebuildable": "reconstruible",
    "kind_expensive": "costoso",
    "kind_results": "resultados",
    "clean_cache":
        "Los cubos: cada espectro de una campaña sobre una misma rejilla de"
        " longitudes de onda, para que una segunda ejecución no los vuelva a"
        " leer todos. Suele ser lo más grande de aquí, y nunca un resultado. Se"
        " reconstruye en minutos por campaña.",
    "clean_spill":
        "Los temporales mapeados de un ajuste demasiado grande para la"
        " memoria. Nada los lee una vez terminado el ajuste, así que lo que"
        " siga aquí pertenece a una ejecución que se interrumpió.",
    "clean_results":
        "Los espectros corregidos, los informes y los ajustes. Aquello para lo"
        " que existe todo esto; nunca se ofrece para borrar.",
    "clean_lbl_science":
        "Los enlaces a través de los cuales mide LBL, uno por exposición. Los"
        " rehace la etapa lbl. Son enlaces simbólicos, así que esto no libera"
        " casi nada y no toca los espectros a los que apuntan.",
    "clean_lbl_plots":
        "Las figuras propias de LBL. Se vuelven a dibujar cada vez que LBL se"
        " ejecuta.",
    "clean_lbl_log": "Los registros de LBL. Solo los lee una persona.",
    "clean_lbl_lblrv":
        "Las tablas de velocidad por línea de LBL, un archivo por exposición,"
        " y lo más grande que escribe LBL. Solo se rehacen volviendo a ejecutar"
        " LBL, lo que toma horas.",
    "clean_lbl_templates":
        "Las plantillas que construyó LBL, una por objeto y ejecución. Solo se"
        " rehacen volviendo a ejecutar la etapa de plantillas de LBL.",
    "clean_lbl_masks":
        "Las máscaras de líneas que construyó LBL, una por objeto y ejecución.",
    "clean_lbl_models": "Los modelos de LBL.",
    "clean_lbl_calib": "Las calibraciones de LBL.",
    "clean_lbl_lblreftable":
        "Las tablas de referencia de LBL, una por objeto y ejecución.",
    "clean_lbl_lblrdb":
        "Las velocidades: el rdb del que sale cada página de VR de cada"
        " informe. Un resultado, y lo más pequeño de esta lista; borrarlo no"
        " libera casi nada y tira la medición.",
    "clean_pycache":
        "Python compilado, que se rehace la próxima vez que se importa el"
        " paquete.",
    "quit": "Salir",
    "quit_title": "salir de pca2d-preclean",
    "quit_yes": "Salir de todos modos",
    "quit_no": "Quedarse",
    "quit_running":
        "Hay una ejecución en curso, y es un subproceso de esta ventana: salir"
        " la detiene. Lo que escribieron las etapas anteriores se queda donde"
        " está; la etapa en curso se pierde. ¿Salir de todos modos?",
    "help_quit_button":
        "Cierra la ventana. Los ajustes se guardan en cada cambio, así que no"
        " se pierde nada al salir; una ejecución en curso se detiene, y la"
        " ventana pregunta antes de hacerlo.",
    "pick_root":
        "elija una raíz de datos: Examinar, junto al campo de arriba",
    "command_pending":
        "marque un objetivo y el comando aparecerá aquí, completo, antes de"
        " ejecutarse",
    "log_pick_root":
        "todavía no hay raíz de datos. Busque la carpeta que contiene UNA"
        " CARPETA POR OBJETIVO de espectros t.fits; nunca se escribe nada en"
        " ella. La raíz de salida se propone a su lado una vez elegida, y la"
        " configuración es la que vino con esta instalación.",
    "log_other_clone":
        "la configuración que usted indica está en OTRA copia de este paquete"
        " (%s), mientras que el código en ejecución es %s. La ejecución usa el"
        " código en ejecución; los ajustes mostrados son los de esa otra copia."
        " Apunte la configuración al mismo lugar, salvo que quiera mezclarlos.",
    "log_no_config":
        "esta instalación no trae config.yaml: busque uno, o abra la ventana"
        " desde un clon del repositorio.",
    "lbl_title": "ajustes de LBL",
    "nothing_export":
        "todos los ajustes son los de la configuración, así que un archivo de"
        " variante no diría nada.",
    "export_title": "guardar estos ajustes como variante",
    "mixed_title": "dos instrumentos",
    "mixed":
        "%s no provienen de un mismo instrumento (%s).\n\nUna ejecución es un"
        " dominio, una rejilla y un conjunto de extensiones, todos leídos del"
        " instrumento, así que objetos de dos espectrógrafos no pueden"
        " ajustarse juntos. Marque solo los de un mismo instrumento.",
    "defaults_title": "guardar como predeterminados",
    "defaults_ask":
        "¿Escribir estos valores en %s como predeterminados para todas las"
        " ejecuciones?\n\n%s\n\nSe conservan los comentarios del archivo.",
    "log_index": "índice de esta raíz de datos: %s",
    "log_out_proposed": "raíz de salida propuesta, junto a los datos: %s",
    "log_scan_start": "leyendo la raíz de datos %s",
    "log_watch_added": "nuevo en la raíz de datos: %s. Leyéndolo.",
    "log_watch_gone": "ya no está en la raíz de datos: %s",
    "log_watch_grew": "más espectros en %s que hace un momento. Leyéndolos.",
    "log_scan_done":
        "%s: %d objetos, %d espectros leídos, %d ya conocidos, %d ausentes",
    "log_scan_none": "ninguna carpeta con espectros bajo %s",
    "log_busy":
        "todavía leyendo la raíz de datos: espere a que termine",
    "log_no_root": "no es una carpeta: %s",
    "log_no_object": "ningún objeto marcado: marque al menos uno en la lista",
    "log_mixed":
        "dos instrumentos marcados (%s): una ejecución es un solo instrumento",
    "log_mixed_refused":
        "se rechaza la ejecución: %s provienen de dos instrumentos (%s)",
    "log_nights": "%s: %d noches en común, de %s",
    "log_nights_thin":
        "%s: %d noches en común de %s. Pocas es una VENTAJA aquí: las"
        " componentes del observador son un parásito siempre presente, lo que"
        " se mide es su patrón, y las noches que ninguna otra estrella vio"
        " amplían el rango de condiciones sobre el que se mide ese patrón.",
    "log_command": "ejecutando: %s",
    "log_ended": "la ejecución terminó, código de salida %d",
    "log_stopping": "pidiendo a la ejecución que se detenga",
    "log_export": "variante escrita: %s",
    "log_defaults": "predeterminados escritos en %s: %s",
    "log_nothing":
        "nada que escribir: todos los ajustes son los de la configuración",
    "log_saved_log": "registro escrito: %s",
    "log_no_outputs": "todavía no hay nada escrito allí: %s",
    "log_no_report":
        "todavía no hay PDF de compilación en %s. Lo escribe la etapa de"
        " figuras, y la etapa LBL le añade las páginas de velocidades",
    "log_report_launched":
        "los ajustes han cambiado desde entonces; este es el informe de la"
        " ejecución que se lanzó: %s",
    "log_failed": "no se pudo iniciar: %s",
    "log_no_command": "no se instaló nada, así que no hay nada que ejecutar",
    "log_installing": "instalando el comando en este entorno, desde %s",
    "no_command_title": "el comando no está instalado",
    "no_command":
        "no se encuentra pca2d-preclean para\n\n  %s\n\n¿Instalarlo ahora"
        " desde %s?\n\nEso ejecuta `pip install -e . --no-deps` allí, lo que"
        " pone el comando en el bin de este intérprete y no toca nada más.",
    "no_command_still": "sigue sin encontrarse tras la instalación",
    "opt_n_star": "componentes estelares",
    "opt_n_earth": "componentes del observador",
    "opt_mean": "parte estática",
    "opt_velocity_term": "ajustar una velocidad por exposición",
    "opt_iters": "barridos como máximo",
    "opt_shrink": "dividir solo lo significativo",
    "opt_weight": "métrica del ajuste de la corrección",
    "help_weight":
        "La métrica en la que se miden las amplitudes de la corrección."
        " `flux`, la nominal: cada muestra tal como la vio el ajuste."
        " `velocity`: cada muestra ponderada por la derivada de la propia"
        " estrella en ese punto, (dT/dv)^2, porque lo que un contaminante le"
        " hace a una velocidad radial es su solapamiento con esa derivada, y un"
        " contaminante plano donde la estrella tiene estructura no mueve"
        " ninguna línea. Implica volver a ajustar las amplitudes de cada"
        " exposición, así que la corrección es más lenta; no cambia nada más:"
        " ni lo que se divide, ni qué muestras se anulan, ni la contracción."
        " Medido en TOI-2120, donde la corrección gana un factor tres, las dos"
        " son indistinguibles: 15.4 +- 1.2 frente a 15.0 +- 1.2 m/s. Los"
        " objetivos donde la corrección CUESTA son los que decidirán.",
    "opt_width_kms": "paso alto (km/s)",
    "opt_dv": "paso de la rejilla (km/s)",
    "opt_nightly_stack": "combinar cada noche",
    "opt_run": "ejecutar LBL (horas)",
    "opt_lbl_prepare": "escribir el árbol de LBL",
    "opt_lbl_before": "medir los espectros entregados",
    "opt_lbl_after": "medir los espectros corregidos",
    "opt_lbl_star_template": "nuestra estrella como plantilla de LBL",
    "opt_lbl_strpca": "componentes estelares extra como RESPROJ",
    "opt_lbl_suffix": "nombre del corregido",
    "opt_lbl_teff": "temperatura efectiva",
    "opt_lbl_template": "archivo de plantilla",
    "opt_lbl_steps": "etapas",
    "opt_lbl_link": "espectros dentro como",
    "opt_lbl_env": "qué LBL",
    "help_lbl_env":
        "Qué LBL mide, por el entorno conda en el que corre su script."
        " `lbl-rapide`, el predeterminado, es la rama rápida de LBL"
        " (test-speed-260918-110104) en un entorno propio, mucho más rápida,"
        " cuyos commits anuncian cada cambio sin efecto en las salidas. `current` es el LBL instalado junto"
        " a este programa, el que usaron todas las ejecuciones anteriores."
        " Una ejecución que no encuentra el elegido se detiene al principio,"
        " antes del ajuste, y dice cómo crearlo.",
    "help_lbl_prepare":
        "Escribir lbl_config.yaml y run_lbl.py junto a las salidas de la"
        " ejecución y poner ambos conjuntos de espectros en las carpetas"
        " science de LBL. Desactivado, la corrección se escribe igual y no se"
        " prepara nada para medir velocidades.",
    "help_lbl_before":
        "Medir también los espectros ENTREGADOS, como su propio objeto en el"
        " mismo árbol LBL. Desactivado, la corrección no se compara con nada y"
        " el resultado no puede leerse como mejor o peor que no hacer nada.",
    "help_lbl_after":
        "Medir los espectros CORREGIDOS. Desactivado con `before` activado, LBL"
        " mide solo los entregados, que es como se construye una vez una serie"
        " de referencia que luego reutilizan todas las ejecuciones.",
    "help_lbl_star_template":
        "Darle a LBL el espectro estelar que construyó este ajuste, en lugar de"
        " dejar que LBL construya el suyo a partir de los espectros corregidos."
        " DESACTIVADO es lo nominal: medido en Proxima, nuestra plantilla hizo"
        " que LBL ajustara líneas un 20% más anchas y duplicó el error por"
        " exposición, de 0.97 a 1.67 m/s, para un rms de 3.34 m/s frente a"
        " 3.04. Lo que le falte a nuestra plantilla, LBL lo paga línea por"
        " línea.",
    "help_lbl_strpca":
        "Con dos o más componentes estelares, las posteriores a la primera se"
        " le dan a LBL como tablas RESPROJ, igual que sus propios gradientes"
        " DTEMP, de modo que el rdb lleva la proyección de cada una por"
        " exposición y se puede buscar una correlación en lugar de suponerla"
        " ausente.",
    "help_lbl_suffix":
        "Cómo se llama el objeto corregido junto al entregado: TOI-2120 y"
        " TOI-2120_PCA2D_2-7. `{tag}` es el número de componentes de la"
        " ejecución, y omitirlo hace que dos ejecuciones escriban sus espectros"
        " corregidos en UNA misma carpeta, donde LBL mide la mezcla sin"
        " avisar.",
    "help_lbl_teff":
        "La temperatura efectiva que se le indica a LBL, que es como elige su"
        " lista de líneas. `auto` la lee del encabezado, y un número la"
        " reemplaza.",
    "help_lbl_template":
        "Un ARCHIVO de plantilla para que LBL lo use, por ruta, en lugar del"
        " que construiría. Vacío es lo nominal: LBL construye una plantilla a"
        " partir de los espectros que mide, que es la serie contra la que se"
        " mide cada conjunto.",
    "help_lbl_steps":
        "Qué etapas de LBL ejecutar, en orden: template, mask, compute,"
        " compile. Menos sirve para retomar un árbol que ya tiene las"
        " anteriores, nunca para saltarse trabajo que necesita una etapa"
        " posterior.",
    "help_lbl_link":
        "Cómo llegan los espectros a las carpetas science de LBL. `symlink`,"
        " el predeterminado, pone allí un enlace: nada se almacena dos veces, y"
        " LBL solo los lee. `copy` pone allí cada espectro por segunda vez,"
        " decenas de gigabytes por campaña, y una carpeta que ya no necesita el"
        " disco de datos. Un disco que no admite enlaces recibe copias diga lo"
        " que diga este ajuste, y la ejecución lo dice al principio.",
    "help_col_snr":
        "La mediana, sobre las exposiciones de este objetivo, de la señal a"
        " ruido de extracción por orden que el pipeline escribió en cada"
        " archivo, para distinguir un objetivo brillante de uno débil antes de"
        " ejecutar nada. Se lee solo de los encabezados y se recuerda, así que"
        " una carpeta se abre al instante la segunda vez.",
    "help_col_exptime":
        "El tiempo de exposición mediano de los archivos de este objetivo, en"
        " segundos. Junto con la S/R y el número de archivos, dice qué tipo de"
        " campaña es: muchas exposiciones cortas de una estrella brillante, o"
        " pocas largas de una débil.",
    "help_col_mag":
        "Haga clic en un encabezado de columna para ordenar por ella, otra vez"
        " para invertir, y dos veces en el encabezado de objeto para volver a"
        " instrumento y luego nombre. Ordenar por ESTA compara bandas que no son"
        " la misma: NIRPS escribe J y SPIRou H, así que una J de 5.3 y una H de"
        " 10.5 se ordenan como números y no como brillos.\n\nEl brillo del"
        " objetivo tal como lo registró SU PROPIO pipeline, con la banda en que"
        " está: NIRPS escribe la magnitud J, SPIRou escribe H, y las dos"
        " difieren en cerca de una magnitud en una enana M, así que se muestra"
        " la banda en lugar de suponerla. Se lee de los encabezados, nunca de"
        " un catálogo.",
    "help_check":
        "Marque un objetivo para ejecutarlo. Varios marcados se ajustan JUNTOS"
        " contra una sola base del observador. Haga clic en la casilla, doble"
        " clic en la fila, o pulse la barra espaciadora. Todos deben venir del"
        " mismo instrumento.",
    "help_all_button":
        "Marca todos los objetivos de la lista, o los desmarca todos.",
    "help_instrument_filter":
        "Muestra u oculta los objetivos de un instrumento. Una ejecución es UN"
        " instrumento, así que ocultar los demás es la vía más rápida hacia una"
        " selección que realmente se puede ejecutar; por la misma razón, las"
        " filas se agrupan y se colorean por instrumento. Ocultar un"
        " instrumento desmarca sus objetivos, y no se borra nada: marque de"
        " nuevo la casilla y vuelven tal como estaban.",
    "help_savedefaults_button":
        "Escribe los ajustes que DIFIEREN de la configuración en el propio"
        " config.yaml, como los predeterminados de los que parten todas las"
        " ejecuciones posteriores. Se conservan los comentarios del archivo, ya"
        " que son las mediciones que eligieron cada valor. Un archivo de"
        " variante deja lo nominal intacto; esto lo cambia.",
    "help_lblwin_button":
        "El bloque LBL en su propia ventana: qué se mide, cómo se llama el"
        " objeto corregido, contra qué plantilla se mide, qué etapas de LBL se"
        " ejecutan. Esa es la etapa que produce velocidades.",
    "help_data_dir":
        "La RAÍZ de entrada, no la carpeta de un objeto: los objetos de abajo"
        " son sus subcarpetas. Nunca se escribe nada en ella.",
    "help_config":
        "config.yaml: todo lo que la ventana no muestra. En él se combinan tres"
        " capas: lo general, lo propio del espectrógrafo y lo propio de un"
        " objetivo.",
    "help_out_dir":
        "Dónde escribe una ejecución. Se propone como `corrected` JUNTO A la"
        " raíz de datos, ya que los espectros corregidos son una copia de la"
        " campaña, decenas de gigabytes, y corresponden al disco donde ya está"
        " la campaña. Vacío usa la raíz de salida de la configuración. Una"
        " ejecución pone su configuración resuelta, su informe, sus espectros"
        " corregidos y sus carpetas LBL bajo <raíz>/<objeto>/<M>-<N>/.",
    "help_rescan":
        "Vuelve a leer la raíz de datos. Lo leído antes se recuerda en un"
        " índice bajo su carpeta personal, nunca en la raíz de datos, así que"
        " solo se leen los archivos añadidos o reemplazados: úselo después de"
        " copiar espectros nuevos.",
    "help_objects":
        "Un objeto marcado: una ejecución individual. VARIOS: se ajustan"
        " juntos contra UNA base del observador, conservando cada uno su propio"
        " espectro estelar por paridad de orden. La atmósfera y el instrumento"
        " son compartidos, las estrellas no, así que una base ajustada sobre"
        " varias estrellas no puede seguir a ninguna de ellas. Deben venir del"
        " mismo instrumento.",
    "help_n_star":
        "Componentes ajustadas en el marco en reposo de la ESTRELLA, además del"
        " propio espectro estelar, que siempre se retira con un coeficiente"
        " exactamente igual a uno. 0 es lo nominal: una amplitud libre delante"
        " de un término en un logaritmo es un exponente sobre el flujo, y un"
        " exponente sobre el espectro medio de una estrella no describe ninguna"
        " estrella. Quitarla llevó TOI-4552 de 14.6 a 11.2 m/s de dispersión"
        " robusta, y TOI-2120 de 14.4 a 14.1.",
    "help_n_earth":
        "Componentes ajustadas en el marco del OBSERVADOR: la atmósfera y el"
        " instrumento. Tres es lo nominal. Más libertad describe mejor el cielo"
        " y se lleva más de la estrella allí donde los dos marcos son"
        " degenerados, que es lo que produce una cobertura de BERV estrecha.",
    "help_mean":
        "La parte del modelo que no tiene amplitud propia: el espectro medio."
        " Hay uno por PARIDAD DE ORDEN, ya que los órdenes pares e impares ven"
        " una longitud de onda con resoluciones distintas, y la cuestión es en"
        " qué marco vive. `star`, lo nominal: un espectro por paridad en el"
        " marco de la ESTRELLA, una mediana por intervalos de BERV retirada una"
        " vez antes de cualquier componente con un coeficiente exactamente"
        " igual a 1, y ninguna media en el marco del observador, de modo que la"
        " corrección divide solo el bloque del observador. `offset`: ningún"
        " espectro estelar, y solo la parte de la media del marco del"
        " observador que DIFIERE entre paridades, dejando la parte compartida"
        " para que la describa el bloque del observador; vuelve a la"
        " corrección. `full`: toda la media del marco del observador por"
        " paridad, incluida la parte compartida. `iterate`: una media por"
        " paridad en CADA marco, reestimada en cada barrido. Medido: TOI-2120,"
        " 20.4 m/s con `star` frente a 31.7 con `offset`; en Proxima la media"
        " del marco del observador por sí sola inyectó 46 m/s y el bloque del"
        " observador por sí solo 48, su suma 19. `iterate` converge con datos"
        " sintéticos y no lo hizo con una campaña completa.",
    "help_velocity_term":
        "Ajusta una velocidad por exposición junto a las componentes, para"
        " mantener el movimiento propio de la estrella fuera del bloque del"
        " observador. Se ajusta y se escribe en los archivos corregidos, nunca"
        " se divide. Medido en TOI-4552 no cambió nada: 16.4 frente a 15.9"
        " m/s.",
    "help_iters":
        "Barridos como máximo. Un barrido resuelve las amplitudes de cada"
        " exposición y luego actualiza cada base contra el residuo de la otra."
        " El ajuste conserva la mejor iteración y se detiene cuando chi2 ha"
        " dejado claramente de bajar, así que esto es un techo, no una"
        " duración.",
    "help_shrink":
        "Divide cada componente del observador solo donde los datos la"
        " detectan, en cada columna, contando todas las exposiciones juntas: un"
        " patrón a media sigma en cada una de N exposiciones se detecta a unos"
        " 0.5 sqrt(N) y se conserva, mientras que una columna donde la"
        " componente es ruido se deja intacta. Sin esto, TOI-4552 perdía 3 m/s"
        " más: 18.1 frente a 15.1.",
    "help_width_kms":
        "El ancho del filtro de Savitzky-Golay que quita el continuo, en km/s,"
        " para que trate una línea igual sea cual sea el paso de la rejilla."
        " Forma parte del cubo: cambiarlo construye uno nuevo.",
    "help_dv":
        "El paso de la rejilla, en km/s. La rejilla es uniforme en logaritmo de"
        " la longitud de onda, lo que convierte un desplazamiento Doppler en"
        " una traslación exacta. Más fino es más pesado: el cubo y el ajuste"
        " crecen con él. Cambiarlo construye un cubo nuevo.",
    "help_nightly_stack":
        "Combina las exposiciones de una noche antes del ajuste. Una decisión"
        " de memoria y nunca de modelado: cada exposición se registra primero"
        " con su propio BERV, así que nada se emborrona, y la corrección sigue"
        " resolviendo las amplitudes de cada exposición contra la base fija."
        " `auto` compara la huella del ajuste con la memoria que permite la"
        " configuración.",
    "help_run":
        "Ejecuta LBL una vez corregidos los espectros. Es la etapa que produce"
        " velocidades, y toma horas. Desactivado, se escribe igualmente todo lo"
        " que LBL necesita y la ejecución dice cómo lanzarlo a mano.",
    "help_stage_cube":
        "Lee cada espectro una vez sobre la rejilla común, con sus pesos."
        " Unos minutos, y queda en caché: una segunda ejecución con los mismos"
        " ajustes lo reutiliza.",
    "help_stage_fit":
        "La descomposición en dos marcos. La larga: imprime un R2 después de"
        " cada barrido.",
    "help_stage_figures":
        "Un PDF: la campaña en el tiempo, los parámetros, y luego una página"
        " por ventana de longitud de onda que muestra cada paso del modelo y lo"
        " que divide la corrección.",
    "help_stage_correct":
        "Escribe los t.fits corregidos, uno por exposición, con el bloque del"
        " observador dividido y las muestras sin peso anuladas.",
    "help_stage_lbl":
        "Le entrega a LBL ambos conjuntos de espectros, los entregados y los"
        " corregidos, como dos objetos en un mismo árbol, para que las"
        " velocidades puedan compararse en lugar de creerse.",
    "help_command":
        "Exactamente lo que ejecutará el botón Ejecutar. Cópielo en una"
        " terminal y hace lo mismo: la ventana es un envoltorio, no una segunda"
        " manera de hacer las cosas.",
    "help_run_button":
        "Inicia el comando de arriba. La salida aparece abajo a medida que se"
        " imprime, y la ejecución se puede detener.",
    "help_stop_button":
        "Pide a la ejecución que se detenga. Lo ya escrito se queda: un cubo,"
        " un ajuste y los archivos corregidos son reutilizables, y una etapa"
        " que no terminó simplemente se vuelve a ejecutar la próxima vez.",
    "help_dry_button":
        "Lo resuelve todo, imprime el plan y no toca nada. La manera honesta"
        " de ver a qué cubo apunta un conjunto de opciones antes de gastar una"
        " hora en él.",
    "help_export_button":
        "Escribe los ajustes que DIFIEREN de la configuración como un archivo"
        " de variante, que es como se hace reproducible una ejecución. Un"
        " archivo que repite lo nominal no dice nada, y se rechaza.",
    "help_savelog_button":
        "Escribe en un archivo lo que la ventana ha mostrado.",
    "help_openout_button":
        "Abre la raíz de salida en el explorador de archivos.",
    "help_openpdf_button":
        "Abre el PDF propio de esta ejecución, en el que se reúne todo: la"
        " estrella, los espectros antes y después, las componentes, las"
        " correlaciones, y las páginas de velocidades una vez que LBL las ha"
        " medido. Es <objeto>_<etiqueta>.pdf en la carpeta de la ejecución,"
        " escrito por la etapa de figuras, así que existe una vez que esa etapa"
        " se ha ejecutado.",
    "help_lang":
        "Elige el idioma de la ventana: haga clic en el que desee. Inglés,"
        " francés, español y portugués.",
    "help_apero":
        "APERO, el pipeline que redujo cada espectro que lee esta ventana"
        " (Cook et al. 2022, PASP 134, 114509). Es quien escribió las"
        " extensiones, la solución en longitud de onda, la señal a ruido por"
        " orden y las velocidades baricéntricas que se usan aquí; este paquete"
        " parte de sus t.fits y nunca vuelve a derivar nada de eso.",
    "help_log":
        "El registro propio de la ejecución, con los colores que le daría una"
        " terminal: verde para el progreso, azul para un número, naranja para"
        " algo omitido, rojo para lo que detiene una ejecución.",
}
