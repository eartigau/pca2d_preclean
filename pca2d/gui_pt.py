"""The window in Portuguese, as written in Portugal (added 2026-09-16).

Every key of gui.EN, with the same placeholders in the same order, in the
European usage the flag on its button stands for: ficheiro, ecrã, registo.
"""

PT = {
    "subtitle": "pré-limpeza em dois referenciais, depois LBL. Escolha uma"
                " raiz de dados, escolha objetos, veja o comando e execute-o.",
    "data_dir": "raiz de dados",
    "config": "configuração",
    "fits_dir": "disco dos produtos (opcional)",
    "help_fits_dir":
        "Onde se GUARDAM os produtos de uma execução: cada pasta de execução"
        " passa a ser uma ligação para este disco, pelo que os espetros"
        " corrigidos, o ajuste e as figuras ficam aí e não no disco interno."
        " Vazio, tudo fica sob a raiz de saída, o que funciona. Um caminho"
        " existe numa ÚNICA máquina, por isso este campo é esvaziado sempre que"
        " o que indica não está presente, e uma execução cujo disco falta"
        " para, em vez de encher discretamente o disco interno.",
    "log_no_disk":
        "a configuração indica %s como o disco onde guardar os produtos, e ele"
        " não está presente: o campo fica vazio, pelo que esta execução"
        " guardaria tudo sob a raiz de saída. Preencha-o se o disco deveria"
        " estar montado.",
    "out_dir": "raiz de saída (opcional)",
    "lbl_dir": "pasta de saída do LBL",
    "help_lbl_dir":
        "A árvore própria do LBL (o seu DATA_DIR): as pastas science que lê, e"
        " os modelos, máscaras, tabelas por linha e rdb que escreve. Proposta"
        " como lbl sob a raiz de saída, e segue essa raiz até que escreva"
        " outra. Uma única árvore para todas as execuções sob uma raiz, já que"
        " o LBL do objeto entregue é o mesmo para todas e demora horas; aponte"
        " para uma árvore existente e essa árvore é usada tal como está.",
    "browse": "Procurar",
    "rescan": "Reler",
    "objects": "objetos",
    "settings": "definições",
    "stages": "etapas",
    "command": "o comando que será executado",
    "output": "saída",
    "col_object": "objeto",
    "col_files": "ficheiros",
    "col_instrument": "instrumento",
    "col_snr": "S/R",
    "col_exptime": "exp. (s)",
    "col_mag": "mag",
    "run_name": "nome da redução",
    "auto": "auto",
    "help_auto_button":
        "Volta a propor um nome a partir dos alvos e das definições tal como"
        " estão: os alvos, e depois seis caracteres de um hash de tudo o que"
        " faz desta redução um resultado diferente. Os mesmos parâmetros dão os"
        " mesmos seis; um número diferente dá outros seis.",
    "timeline": "quando foram observadas as campanhas assinaladas",
    "timeline_none": "assinale um alvo para ver quando foi observado",
    "timeline_note": "a manter de %s a %s: %d exposições de %d",
    "whole": "tudo",
    "min_rjd": "de",
    "max_rjd": "a",
    "exists": "esta execução já existe",
    "help_run_name":
        "Um nome para esta execução. Os seus produtos vão para"
        " <raiz de saída>/_NOME/ e o seu objeto LBL é"
        " <objeto>_PCA2D_<M-N>_NOME, pelo que duas execuções dos mesmos alvos"
        " com definições diferentes nunca escrevem numa mesma pasta nem sob um"
        " mesmo nome LBL, onde o LBL mediria a mistura sem dizer nada. Vazio é"
        " o caminho nominal. É proposto a partir das datas quando estão"
        " definidas, e pode ser qualquer coisa.",
    "help_dates":
        "Mantém apenas as exposições entre estas duas datas, em data juliana"
        " reduzida (BJD - 2400000). O que fica de fora não é ajustado nem"
        " corrigido. É assim que se reduz uma campanha a uma temporada, o que"
        " a cobertura baricêntrica por vezes exige: dois cortes de 14 noites de"
        " TOI-4552 com a mesma relação sinal-ruído diferem por um fator de"
        " sessenta em cobertura, e a correção muda de sinal entre eles.",
    "berv": "cobertura baricêntrica dos alvos assinalados",
    "berv_none": "assinale um alvo para ver a sua cobertura baricêntrica",
    "berv_wait": "ainda não foi lida nenhuma velocidade baricêntrica",
    "berv_building": "EM CONSTRUÇÃO: ainda a ler",
    "berv_note":
        "cobertura efetiva %.0f km/s, extensão %.0f, no máximo %s possível para"
        " a latitude eclíptica deste alvo; intervalos de %.0f km/s,"
        " %d exposições",
    "help_berv":
        "Que parte do intervalo baricêntrico as campanhas assinaladas cobrem"
        " de facto, em intervalos de 3 km/s, uma cor por instrumento e"
        " empilhados. É a grandeza que decide se a correção ajuda: o bloco do"
        " observador só é identificável porque a estrela se desloca através"
        " dele. Dois cortes de 14 noites de TOI-4552 com a mesma relação"
        " sinal-ruído, um com 42.8 km/s de extensão e o outro com 0.7,"
        " passaram de um ganho de um fator dois para uma perda de um fator"
        " 2.4. A cobertura EFETIVA conta apenas os intervalos que contêm uma"
        " exposição, pelo que a uma campanha observada em dois extremos e em"
        " nenhum ponto intermédio não é creditado o intervalo entre eles. O"
        " terceiro número é o que o céu PERMITE: o |BERV| de um alvo nunca"
        " excede 29.78 cos(latitude eclíptica) km/s, pelo que TOI-1452, a"
        " +80.5 graus, só pode abranger 9.8 km/s e nenhuma quantidade de"
        " observações separará os dois referenciais para ele. Uma cobertura"
        " inferior à extensão indica lacunas, que mais noites podem preencher;"
        " uma extensão inferior ao possível indica uma campanha jovem; um"
        " possível pequeno indica que o alvo não é o adequado.",
    "run": "Executar",
    "stop": "Parar",
    "dry": "Simulação",
    "export": "Exportar YAML...",
    "savelog": "Guardar registo...",
    "openout": "Abrir saídas",
    "openpdf": "Abrir PDF de compilação",
    "savedefaults": "Guardar como predefinições...",
    "all": "todos",
    "none": "nenhum",
    "idle": "em espera",
    "running": "em curso",
    "lang": "Português",
    "snr_berv": "relação sinal-ruído em função da velocidade baricêntrica",
    "snr_none": "assinale um alvo para ver onde estão as suas melhores noites",
    "help_snr_berv":
        "Um ponto por espetro, na cor da sua estrela: a relação sinal-ruído que"
        " o APERO mediu em função da velocidade baricêntrica em que foi obtido."
        " O histograma acima diz que velocidades uma campanha cobre; este diz"
        " COM QUE as cobre. O ajuste pondera um espetro por 1/sigma^2, pelo que"
        " um intervalo coberto numa ponta pelas piores noites de uma campanha"
        " não é o intervalo que o ajuste realmente vê, e os dois referenciais"
        " separam-se pior do que a cobertura promete.",
    "tab_targets": "  alvos  ",
    "tab_settings": "  definições  ",
    "tab_lbl": "  LBL  ",
    "tab_run": "  análise  ",
    "tab_clean": "  limpeza  ",
    "tab_runs": "  passagens  ",
    "runs_title": "todas as passagens que esta pasta de saída contém",
    "help_runs":
        "O que já foi tentado, lido do que cada passagem deixou escrito: uma"
        " passagem regista os ajustes que resolveu, o comando que recebeu e a"
        " hora a que começou, por isso esta lista são as próprias passagens e"
        " não um histórico guardado pela janela. Escolha uma para ver todas as"
        " suas opções e abrir o seu PDF de compilação.",
    "runs_rescan": "Reler",
    "help_runs_rescan":
        "Percorrer de novo a pasta de saída e listar todas as passagens que lá"
        " estão. Uma passagem a decorrer aparece assim que tiver escrito a sua"
        " configuração resolvida, ou seja, antes da sua primeira etapa.",
    "runs_open_pdf": "Abrir o seu PDF de compilação",
    "help_runs_open_pdf":
        "Abrir o único PDF da passagem escolhida na lista: as suas figuras, os"
        " seus números e, depois de o LBL medir, as suas velocidades. Escrito"
        " pela etapa das figuras e outra vez depois do LBL.",
    "runs_open_folder": "Abrir a sua pasta",
    "help_runs_open_folder":
        "Mostrar a pasta da passagem escolhida no explorador de ficheiros: a"
        " configuração resolvida, o ajuste, os espectros corrigidos e o"
        " relatório estão lá todos.",
    "runs_scanning": "a ler a pasta de saída...",
    "runs_found": "%d passagem(ns) em %s",
    "runs_none": "ainda não há passagens em %s: uma passagem escreve lá a sua"
                 " configuração resolvida antes da sua primeira etapa",
    "runs_pick": "escolha primeiro uma passagem na lista",
    "runs_no_pdf": "ainda não há PDF de compilação em %s",
    "runs_command": "o comando que recebeu",
    "col_hash": "hash", "col_started": "começou",
    "col_targets": "alvos", "col_tag": "contadores",
    "col_components": "estrela + observador", "col_report": "PDF",
    "col_setting": "ajuste", "col_value": "esta passagem",
    "clean_title": "o que o pca2d deixou nestes discos",
    "help_clean":
        "Cada sítio onde este programa escreve bytes, medido. A cache de"
        " cubos, os temporários de um ajuste demasiado grande para a memória e"
        " os registos do LBL podem ser apagados a qualquer momento: refazem-se"
        " a partir do que continua aqui, e um botão esvazia-os todos. Tudo o"
        " resto também pode ser apagado, linha a linha: escolha-a na lista e"
        " Apagar a seleção. Os ajustes, os modelos, as máscaras e as tabelas"
        " por linha só voltam executando as horas que os produziram, e os"
        " espetros corrigidos, os relatórios e as velocidades só voltando a"
        " executar tudo a partir dos espetros; a janela diz quais dessas horas"
        " está prestes a gastar antes de apagar seja o que for. Uma pasta"
        " movida para outro disco é contada onde realmente está; uma ligação"
        " para os espetros de outra pessoa não pesa nada, porque apagá-la não"
        " liberta nada.",
    "clean_measure": "Medir",
    "clean_purge": "Apagar o que pode sair",
    "clean_selected": "Apagar a seleção",
    "help_clean_measure":
        "Percorre cada pasta e soma o que contém. Num disco de rede isto demora"
        " alguns segundos e corre à parte, pelo que a janela continua"
        " utilizável enquanto conta.",
    "help_clean_purge":
        "Esvazia as pastas marcadas como temporárias ou reconstruíveis, e nada"
        " mais: as que só custam tempo a recuperar. Pergunta antes, e indica o"
        " total que vai libertar. As pastas em si ficam: a próxima execução"
        " espera encontrá-las. Para tudo o resto, escolha as linhas e use"
        " Apagar a seleção.",
    "help_clean_selected":
        "Apaga as linhas escolhidas na lista, sejam de que tipo forem: uma"
        " cache de cubos, os modelos que levaram uma tarde, as próprias"
        " velocidades. Shift ou command escolhe várias. Indica cada pasta e diz"
        " o que custaria recuperá-la antes de apagar seja o que for, e isso é"
        " a única coisa entre si e um disco com espaço livre.",
    "clean_measuring": "a medir…",
    "clean_totals": "%s no total, dos quais %s podem ser libertados",
    "clean_none": "nada medido ainda: carregue em Medir",
    "clean_confirm_title": "apagar os ficheiros reconstruíveis",
    "clean_confirm":
        "Prestes a libertar %s em %d sítios:\n\n%s\n\nNada disto é um"
        " resultado: os cubos voltam a ser lidos a partir dos espetros, e o LBL"
        " volta a escrever os seus registos. Continuar?",
    "clean_go": "Varrer",
    "clean_keep": "Manter",
    "clean_freed": "%s libertados em %d sítios",
    "clean_nothing": "nada a libertar: não há aqui temporários nem cache",
    "clean_pick":
        "escolha primeiro uma linha ou várias na lista: o botão apaga o que"
        " está escolhido",
    "clean_while_running":
        "HÁ UMA EXECUÇÃO EM CURSO, e lê destas pastas. Reconstrói o que"
        " encontrar em falta, pelo que apagar agora custa a essa execução o"
        " tempo de o refazer, a meio da etapa em que está.",
    "clean_confirm_pick_title": "apagar a seleção",
    "clean_confirm_pick":
        "Prestes a apagar %s em %d sítios:\n\n%s\n\n%s\n\nContinuar?",
    "clean_cost_rebuildable":
        "Tudo se refaz a partir do que fica aqui: minutos, e uma execução que"
        " volta a ler os espetros.",
    "clean_cost_expensive":
        "Uma parte só volta executando o que a produziu: um ajuste, ou uma"
        " passagem do LBL por cada exposição. Horas, não minutos.",
    "clean_cost_results":
        "UMA PARTE É UM RESULTADO: espetros corrigidos, um relatório, ou as"
        " velocidades que um rdb contém. Nada aqui os refaz. Só voltar a"
        " executar tudo, a partir dos espetros, e isso são as horas que levou"
        " da primeira vez.",
    "col_size": "tamanho",
    "col_nfiles": "ficheiros",
    "col_kind": "tipo",
    "kind_scratch": "temporário",
    "kind_rebuildable": "reconstruível",
    "kind_expensive": "dispendioso",
    "kind_results": "resultados",
    "clean_cache":
        "Os cubos: cada espetro de uma campanha numa mesma grelha de"
        " comprimentos de onda, para que uma segunda execução não os volte a"
        " ler todos. Normalmente o maior de tudo aqui, e nunca um resultado."
        " Reconstrói-se em minutos por campanha.",
    "clean_spill":
        "Os temporários mapeados de um ajuste demasiado grande para caber em"
        " memória. Nada os lê depois de o ajuste terminar, pelo que o que"
        " ainda estiver aqui pertence a uma execução que foi interrompida.",
    "clean_results":
        "Os espetros corrigidos, os relatórios e os ajustes. Aquilo para que"
        " tudo isto existe; nunca é proposto para apagar.",
    "clean_lbl_science":
        "As ligações através das quais o LBL mede, uma por exposição. Refeitas"
        " pela etapa lbl. São ligações simbólicas, pelo que isto não liberta"
        " quase nada e não toca nos espetros para que apontam.",
    "clean_lbl_plots":
        "As figuras do próprio LBL. Voltam a ser desenhadas sempre que o LBL"
        " corre.",
    "clean_lbl_log": "Os registos do LBL. Só uma pessoa os lê.",
    "clean_lbl_lblrv":
        "As tabelas de velocidade por linha do LBL, um ficheiro por exposição,"
        " e o maior de tudo o que o LBL escreve. Só se refazem voltando a"
        " executar o LBL, o que leva horas.",
    "clean_lbl_templates":
        "Os modelos que o LBL construiu, um por objeto e execução. Só se"
        " refazem voltando a executar a etapa de modelos do LBL.",
    "clean_lbl_masks":
        "As máscaras de linhas que o LBL construiu, uma por objeto e execução.",
    "clean_lbl_models": "Os modelos estelares do LBL.",
    "clean_lbl_calib": "As calibrações do LBL.",
    "clean_lbl_lblreftable":
        "As tabelas de referência do LBL, uma por objeto e execução.",
    "clean_lbl_lblrdb":
        "As velocidades: o rdb de onde sai cada página de VR de cada"
        " relatório. Um resultado, e o mais pequeno desta lista; apagá-lo não"
        " liberta quase nada e deita fora a medição.",
    "clean_pycache":
        "Python compilado, refeito da próxima vez que o pacote é importado.",
    "quit": "Sair",
    "quit_title": "sair do pca2d-preclean",
    "quit_yes": "Sair mesmo assim",
    "quit_no": "Ficar",
    "quit_running":
        "Há uma execução em curso, e é um subprocesso desta janela: sair"
        " para-a. O que as etapas anteriores escreveram fica onde está; a etapa"
        " em curso perde-se. Sair mesmo assim?",
    "help_quit_button":
        "Fecha a janela. As definições são guardadas a cada alteração, pelo"
        " que nada se perde ao sair; uma execução em curso é parada, e a janela"
        " pergunta antes de o fazer.",
    "pick_root":
        "escolha uma raiz de dados: Procurar, ao lado do campo no topo",
    "command_pending":
        "assinale um alvo e o comando aparece aqui, completo, antes de ser"
        " executado",
    "log_pick_root":
        "ainda não há raiz de dados. Procure a pasta que contém UMA PASTA POR"
        " ALVO de espetros t.fits; nunca se escreve nada nela. A raiz de saída"
        " é proposta ao lado dela assim que for escolhida, e a configuração é a"
        " que veio com esta instalação.",
    "log_other_clone":
        "a configuração que indica está NOUTRA cópia deste pacote (%s),"
        " enquanto o código em execução é %s. A execução usa o código em"
        " execução; as definições mostradas são as dessa outra cópia. Aponte a"
        " configuração para o mesmo sítio, a menos que queira misturá-los.",
    "log_no_config":
        "esta instalação não trouxe config.yaml: procure um, ou abra a janela"
        " a partir de um clone do repositório.",
    "lbl_title": "definições do LBL",
    "nothing_export":
        "todas as definições são as da configuração, pelo que um ficheiro de"
        " variante não diria nada.",
    "export_title": "guardar estas definições como variante",
    "mixed_title": "dois instrumentos",
    "mixed":
        "%s não vêm de um mesmo instrumento (%s).\n\nUma execução é um"
        " domínio, uma grelha e um conjunto de extensões, todos lidos a partir"
        " do instrumento, pelo que objetos de dois espetrógrafos não podem ser"
        " ajustados em conjunto. Assinale apenas os de um mesmo instrumento.",
    "defaults_title": "guardar como predefinições",
    "defaults_ask":
        "Escrever estes valores em %s como predefinições para todas as"
        " execuções?\n\n%s\n\nOs comentários do ficheiro são mantidos.",
    "log_index": "índice desta raiz de dados: %s",
    "log_out_proposed": "raiz de saída proposta, ao lado dos dados: %s",
    "log_scan_start": "a ler a raiz de dados %s",
    "log_watch_added": "novo na raiz de dados: %s. A ler.",
    "log_watch_gone": "já não está na raiz de dados: %s",
    "log_watch_grew": "mais espetros em %s do que há pouco. A lê-los.",
    "log_scan_done":
        "%s: %d objetos, %d espetros lidos, %d já conhecidos, %d desaparecidos",
    "log_scan_none": "nenhuma pasta com espetros sob %s",
    "log_busy": "ainda a ler a raiz de dados: aguarde que termine",
    "log_no_root": "não é uma pasta: %s",
    "log_no_object":
        "nenhum objeto assinalado: assinale pelo menos um na lista",
    "log_mixed":
        "dois instrumentos assinalados (%s): uma execução é um só instrumento",
    "log_mixed_refused":
        "execução recusada: %s vêm de dois instrumentos (%s)",
    "log_nights": "%s: %d noites em comum, de %s",
    "log_nights_thin":
        "%s: %d noites em comum de %s. Poucas é uma VANTAGEM aqui: as"
        " componentes do observador são um parasita sempre presente, o que se"
        " mede é o seu padrão, e as noites que nenhuma outra estrela viu"
        " alargam o leque de condições sobre o qual esse padrão é medido.",
    "log_command": "a executar: %s",
    "log_ended": "a execução terminou, código de saída %d",
    "log_stopping": "a pedir à execução que pare",
    "log_export": "variante escrita: %s",
    "log_defaults": "predefinições escritas em %s: %s",
    "log_nothing":
        "nada a escrever: todas as definições são as da configuração",
    "log_saved_log": "registo escrito: %s",
    "log_no_outputs": "ainda nada foi escrito aí: %s",
    "log_no_report":
        "ainda não há PDF de compilação em %s. A etapa de figuras escreve-o, e"
        " a etapa LBL acrescenta-lhe as páginas de velocidades",
    "log_report_launched":
        "as definições mudaram entretanto; este é o relatório da execução que"
        " foi lançada: %s",
    "log_failed": "não foi possível iniciar: %s",
    "log_no_command":
        "nada foi instalado, pelo que não há nada para executar",
    "log_installing": "a instalar o comando neste ambiente, a partir de %s",
    "no_command_title": "o comando não está instalado",
    "no_command":
        "não se encontra o pca2d-preclean para\n\n  %s\n\nInstalá-lo agora a"
        " partir de %s?\n\nIsso executa `pip install -e . --no-deps` nesse"
        " sítio, o que coloca o comando no bin deste interpretador e não toca"
        " em mais nada.",
    "no_command_still": "continua sem ser encontrado após a instalação",
    "opt_n_star": "componentes estelares",
    "opt_n_earth": "componentes do observador",
    "opt_mean": "parte estática",
    "opt_velocity_term": "ajustar uma velocidade por exposição",
    "opt_iters": "varrimentos no máximo",
    "opt_shrink": "dividir só o que é significativo",
    "opt_weight": "métrica do ajuste da correção",
    "help_weight":
        "A métrica em que as amplitudes da correção são medidas. `flux`, a"
        " nominal: cada amostra tal como o ajuste a viu. `velocity`: cada"
        " amostra ponderada pela derivada da própria estrela nesse ponto,"
        " (dT/dv)^2, porque o que um contaminante faz a uma velocidade radial"
        " é a sua sobreposição com essa derivada, e um contaminante plano onde"
        " a estrela tem estrutura não desloca nenhuma linha. Implica voltar a"
        " ajustar as amplitudes de cada exposição, pelo que a correção é mais"
        " lenta; não muda mais nada: nem o que é dividido, nem que amostras são"
        " anuladas, nem a contração. Medido em TOI-2120, onde a correção ganha"
        " um fator três, as duas são indistinguíveis: 15.4 +- 1.2 contra"
        " 15.0 +- 1.2 m/s. Os alvos onde a correção CUSTA são os que vão"
        " decidir.",
    "opt_width_kms": "passa-alto (km/s)",
    "opt_dv": "passo da grelha (km/s)",
    "opt_nightly_stack": "combinar cada noite",
    "opt_run": "executar o LBL (horas)",
    "opt_lbl_prepare": "escrever a árvore do LBL",
    "opt_lbl_before": "medir os espetros entregues",
    "opt_lbl_after": "medir os espetros corrigidos",
    "opt_lbl_star_template": "a nossa estrela como modelo do LBL",
    "opt_lbl_strpca": "componentes estelares extra como RESPROJ",
    "opt_lbl_suffix": "nome do corrigido",
    "opt_lbl_teff": "temperatura efetiva",
    "opt_lbl_template": "ficheiro de modelo",
    "opt_lbl_steps": "etapas",
    "opt_lbl_link": "espetros lá dentro como",
    "help_lbl_prepare":
        "Escrever lbl_config.yaml e run_lbl.py ao lado das saídas da execução"
        " e colocar os dois conjuntos de espetros nas pastas science do LBL."
        " Desligado, a correção é escrita na mesma e nada é preparado para"
        " medir velocidades.",
    "help_lbl_before":
        "Medir também os espetros ENTREGUES, como objeto próprio na mesma"
        " árvore LBL. Desligado, a correção não é comparada com nada e o"
        " resultado não pode ser lido como melhor ou pior do que não fazer"
        " nada.",
    "help_lbl_after":
        "Medir os espetros CORRIGIDOS. Desligado com `before` ligado, o LBL"
        " mede apenas os entregues, que é como se constrói uma vez uma série de"
        " referência que depois todas as execuções reutilizam.",
    "help_lbl_star_template":
        "Dar ao LBL o espetro estelar que este ajuste construiu, em vez de"
        " deixar o LBL construir o seu a partir dos espetros corrigidos."
        " DESLIGADO é o nominal: medido em Proxima, o nosso modelo fez o LBL"
        " ajustar linhas 20% mais largas e duplicou o erro por exposição, de"
        " 0.97 para 1.67 m/s, para um rms de 3.34 m/s contra 3.04. O que falta"
        " ao nosso modelo, o LBL paga linha a linha.",
    "help_lbl_strpca":
        "Com duas ou mais componentes estelares, as posteriores à primeira são"
        " entregues ao LBL como tabelas RESPROJ, tal como os seus próprios"
        " gradientes DTEMP, pelo que o rdb traz a projeção de cada uma por"
        " exposição e pode procurar-se uma correlação em vez de a supor"
        " ausente.",
    "help_lbl_suffix":
        "Como se chama o objeto corrigido ao lado do entregue: TOI-2120 e"
        " TOI-2120_PCA2D_2-7. `{tag}` é o número de componentes da execução, e"
        " omiti-lo faz com que duas execuções escrevam os seus espetros"
        " corrigidos numa MESMA pasta, onde o LBL mede a mistura sem avisar.",
    "help_lbl_teff":
        "A temperatura efetiva indicada ao LBL, que é como ele escolhe a sua"
        " lista de linhas. `auto` lê-a do cabeçalho, e um número substitui-a.",
    "help_lbl_template":
        "Um FICHEIRO de modelo para o LBL usar, por caminho, em vez do que ele"
        " construiria. Vazio é o nominal: o LBL constrói um modelo a partir dos"
        " espetros que mede, que é a série contra a qual cada conjunto é"
        " medido.",
    "help_lbl_steps":
        "Que etapas do LBL executar, por ordem: template, mask, compute,"
        " compile. Menos serve para retomar uma árvore que já tem as"
        " anteriores, nunca para saltar trabalho de que uma etapa posterior"
        " precisa.",
    "help_lbl_link":
        "Como os espetros chegam às pastas science do LBL. `symlink`, a"
        " predefinição, coloca lá uma ligação: nada é guardado duas vezes, e o"
        " LBL só os lê. `copy` coloca lá cada espetro uma segunda vez, dezenas"
        " de gigabytes por campanha, e uma pasta que já não precisa do disco de"
        " dados. Um disco que não aceita ligações recebe cópias seja qual for"
        " esta definição, e a execução di-lo logo no início.",
    "help_col_snr":
        "A mediana, sobre as exposições deste alvo, da relação sinal-ruído de"
        " extração por ordem que o pipeline escreveu em cada ficheiro, para"
        " distinguir um alvo brilhante de um fraco antes de executar seja o"
        " que for. Lida apenas dos cabeçalhos e memorizada, pelo que uma pasta"
        " abre de imediato da segunda vez.",
    "help_col_exptime":
        "O tempo de exposição mediano dos ficheiros deste alvo, em segundos."
        " Com a S/R e o número de ficheiros, diz que tipo de campanha é esta:"
        " muitas exposições curtas de uma estrela brilhante, ou poucas longas"
        " de uma fraca.",
    "help_col_mag":
        "Clique num cabeçalho de coluna para ordenar por ela, outra vez para"
        " inverter, e duas vezes no cabeçalho de objeto para voltar a"
        " instrumento e depois nome. Ordenar por ESTA compara bandas que não"
        " são a mesma: o NIRPS escreve J e o SPIRou H, pelo que um J de 5.3 e"
        " um H de 10.5 são ordenados como números e não como brilhos.\n\nO"
        " brilho do alvo tal como o registou o SEU PRÓPRIO pipeline, com a"
        " banda em que está: o NIRPS escreve a magnitude J, o SPIRou escreve H,"
        " e as duas diferem de cerca de uma magnitude numa anã M, pelo que a"
        " banda é mostrada em vez de suposta. Lido dos cabeçalhos, nunca de um"
        " catálogo.",
    "help_check":
        "Assinale um alvo para o executar. Vários assinalados são ajustados EM"
        " CONJUNTO contra uma única base do observador. Clique na caixa, faça"
        " duplo clique na linha, ou carregue na barra de espaços. Têm de vir"
        " todos do mesmo instrumento.",
    "help_all_button":
        "Assinala todos os alvos da lista, ou retira a marca a todos.",
    "help_instrument_filter":
        "Mostra ou oculta os alvos de um instrumento. Uma execução é UM"
        " instrumento, pelo que ocultar os outros é o caminho mais rápido para"
        " uma seleção que se pode realmente executar; pela mesma razão, as"
        " linhas são agrupadas e coloridas por instrumento. Ocultar um"
        " instrumento retira a marca aos seus alvos, e nada é apagado: volte a"
        " assinalar a caixa e regressam como estavam.",
    "help_savedefaults_button":
        "Escreve as definições que DIFEREM da configuração no próprio"
        " config.yaml, como as predefinições de que partem todas as execuções"
        " seguintes. Os comentários do ficheiro são mantidos, já que são as"
        " medições que escolheram cada valor. Um ficheiro de variante deixa o"
        " nominal intacto; isto altera-o.",
    "help_lblwin_button":
        "O bloco LBL na sua própria janela: o que é medido, como se chama o"
        " objeto corrigido, contra que modelo é medido, que etapas do LBL são"
        " executadas. É essa a etapa que produz velocidades.",
    "help_data_dir":
        "A RAIZ de entrada, não a pasta de um objeto: os objetos abaixo são as"
        " suas subpastas. Nunca se escreve nada nela.",
    "help_config":
        "config.yaml: tudo o que a janela não mostra. Nele combinam-se três"
        " camadas: o que é geral, o que é próprio do espetrógrafo e o que é"
        " próprio de um alvo.",
    "help_out_dir":
        "Onde uma execução escreve. Proposta como `corrected` AO LADO da raiz"
        " de dados, já que os espetros corrigidos são uma cópia da campanha,"
        " dezenas de gigabytes, e pertencem ao disco onde a campanha já está."
        " Vazio usa a raiz de saída da configuração. Uma execução coloca a sua"
        " configuração resolvida, o seu relatório, os seus espetros corrigidos"
        " e as suas pastas LBL sob <raiz>/<objeto>/<M>-<N>/.",
    "help_rescan":
        "Volta a ler a raiz de dados. O que foi lido antes é memorizado num"
        " índice na sua pasta pessoal, nunca na raiz de dados, pelo que só são"
        " lidos os ficheiros acrescentados ou substituídos: use-o depois de"
        " copiar espetros novos.",
    "help_objects":
        "Um objeto assinalado: uma execução individual. VÁRIOS: são ajustados"
        " em conjunto contra UMA base do observador, mantendo cada um o seu"
        " próprio espetro estelar por paridade de ordem. A atmosfera e o"
        " instrumento são partilhados, as estrelas não, pelo que uma base"
        " ajustada sobre várias estrelas não pode seguir nenhuma delas. Têm de"
        " vir do mesmo instrumento.",
    "help_n_star":
        "Componentes ajustadas no referencial de repouso da ESTRELA, além do"
        " próprio espetro estelar, que é sempre retirado com um coeficiente"
        " exatamente igual a um. 0 é o nominal: uma amplitude livre à frente de"
        " um termo num logaritmo é um expoente sobre o fluxo, e um expoente"
        " sobre o espetro médio de uma estrela não descreve nenhuma estrela."
        " Retirá-la levou TOI-4552 de 14.6 para 11.2 m/s de dispersão robusta,"
        " e TOI-2120 de 14.4 para 14.1.",
    "help_n_earth":
        "Componentes ajustadas no referencial do OBSERVADOR: a atmosfera e o"
        " instrumento. Três é o nominal. Mais liberdade descreve melhor o céu e"
        " leva consigo mais da estrela onde os dois referenciais são"
        " degenerados, que é o que uma cobertura de BERV estreita provoca.",
    "help_mean":
        "A parte do modelo que não tem amplitude própria: o espetro médio. Há"
        " um por PARIDADE DE ORDEM, já que as ordens pares e ímpares veem um"
        " comprimento de onda com resoluções diferentes, e a questão é em que"
        " referencial vive. `star`, o nominal: um espetro por paridade no"
        " referencial da ESTRELA, uma mediana por intervalos de BERV retirada"
        " uma vez antes de qualquer componente com um coeficiente exatamente"
        " igual a 1, e nenhuma média no referencial do observador, pelo que a"
        " correção divide apenas o bloco do observador. `offset`: nenhum"
        " espetro estelar, e apenas a parte da média no referencial do"
        " observador que DIFERE entre paridades, ficando a parte partilhada"
        " para ser descrita pelo bloco do observador; regressa à correção."
        " `full`: toda a média no referencial do observador por paridade,"
        " incluindo a parte partilhada. `iterate`: uma média por paridade em"
        " CADA referencial, reestimada a cada varrimento. Medido: TOI-2120,"
        " 20.4 m/s com `star` contra 31.7 com `offset`; em Proxima a média no"
        " referencial do observador sozinha injetou 46 m/s e o bloco do"
        " observador sozinho 48, a sua soma 19. `iterate` converge em dados"
        " sintéticos e não convergiu numa campanha completa.",
    "help_velocity_term":
        "Ajusta uma velocidade por exposição ao lado das componentes, para"
        " manter o movimento próprio da estrela fora do bloco do observador. É"
        " ajustada e escrita nos ficheiros corrigidos, nunca dividida. Medido"
        " em TOI-4552 não mudou nada: 16.4 contra 15.9 m/s.",
    "help_iters":
        "Varrimentos no máximo. Um varrimento resolve as amplitudes de cada"
        " exposição e depois atualiza cada base contra o resíduo da outra. O"
        " ajuste mantém a melhor iteração e para quando o chi2 claramente"
        " deixou de descer, pelo que isto é um teto, não uma duração.",
    "help_shrink":
        "Divide cada componente do observador apenas onde os dados a detetam,"
        " em cada coluna, contando todas as exposições em conjunto: um padrão"
        " a meio sigma em cada uma de N exposições é detetado a cerca de"
        " 0.5 sqrt(N) e é mantido, enquanto uma coluna onde a componente é"
        " ruído fica intacta. Sem isto, TOI-4552 perdia mais 3 m/s: 18.1"
        " contra 15.1.",
    "help_width_kms":
        "A largura do filtro de Savitzky-Golay que retira o contínuo, em km/s,"
        " para que trate uma linha da mesma forma seja qual for o passo da"
        " grelha. Faz parte do cubo: alterá-la constrói um novo.",
    "help_dv":
        "O passo da grelha, em km/s. A grelha é uniforme no logaritmo do"
        " comprimento de onda, o que torna um desvio Doppler uma translação"
        " exata. Mais fino é mais pesado: o cubo e o ajuste crescem com ele."
        " Alterá-lo constrói um cubo novo.",
    "help_nightly_stack":
        "Combina as exposições de uma noite antes do ajuste. Uma decisão de"
        " memória e nunca de modelação: cada exposição é primeiro registada com"
        " o seu próprio BERV, pelo que nada fica esborratado, e a correção"
        " continua a resolver as amplitudes de cada exposição contra a base"
        " fixa. `auto` compara a pegada do ajuste com a memória que a"
        " configuração permite.",
    "help_run":
        "Executa o LBL assim que os espetros estão corrigidos. É a etapa que"
        " produz velocidades, e leva horas. Desligado, tudo o que o LBL precisa"
        " é escrito na mesma e a execução diz como o lançar à mão.",
    "help_stage_cube":
        "Lê cada espetro uma vez para a grelha comum, com os seus pesos."
        " Alguns minutos, e fica em cache: uma segunda execução com as mesmas"
        " definições reutiliza-o.",
    "help_stage_fit":
        "A decomposição em dois referenciais. A longa: imprime um R2 após cada"
        " varrimento.",
    "help_stage_figures":
        "Um PDF: a campanha ao longo do tempo, os parâmetros, e depois uma"
        " página por janela de comprimento de onda que mostra cada passo do"
        " modelo e o que a correção divide.",
    "help_stage_correct":
        "Escreve os t.fits corrigidos, um por exposição, com o bloco do"
        " observador dividido e as amostras sem peso anuladas.",
    "help_stage_lbl":
        "Entrega ao LBL os dois conjuntos de espetros, os entregues e os"
        " corrigidos, como dois objetos numa mesma árvore, para que as"
        " velocidades possam ser comparadas em vez de acreditadas.",
    "help_command":
        "Exatamente o que o botão Executar vai correr. Copie-o para um"
        " terminal e faz o mesmo: a janela é um invólucro, não uma segunda"
        " maneira de fazer as coisas.",
    "help_run_button":
        "Inicia o comando acima. A saída aparece abaixo à medida que é"
        " impressa, e a execução pode ser parada.",
    "help_stop_button":
        "Pede à execução que pare. O que já foi escrito fica: um cubo, um"
        " ajuste e os ficheiros corrigidos são reutilizáveis, e uma etapa que"
        " não terminou simplesmente corre de novo da próxima vez.",
    "help_dry_button":
        "Resolve tudo, imprime o plano e não toca em nada. A forma honesta de"
        " ver para que cubo aponta um conjunto de opções antes de gastar uma"
        " hora com ele.",
    "help_export_button":
        "Escreve as definições que DIFEREM da configuração como um ficheiro de"
        " variante, que é como uma execução se torna reprodutível. Um ficheiro"
        " que repete o nominal não diz nada, e é recusado.",
    "help_savelog_button":
        "Escreve num ficheiro o que a janela mostrou.",
    "help_openout_button":
        "Abre a raiz de saída no explorador de ficheiros.",
    "help_openpdf_button":
        "Abre o PDF próprio desta execução, onde tudo é reunido: a estrela, os"
        " espetros antes e depois, as componentes, as correlações, e as páginas"
        " de velocidades assim que o LBL as tiver medido. É"
        " <objeto>_<etiqueta>.pdf na pasta da execução, escrito pela etapa de"
        " figuras, pelo que existe assim que essa etapa tiver corrido.",
    "help_lang":
        "Escolhe a língua da janela: clique na que pretende. Inglês, francês,"
        " espanhol e português.",
    "help_apero":
        "O APERO, o pipeline que reduziu cada espetro que esta janela lê (Cook"
        " et al. 2022, PASP 134, 114509). Foi ele que escreveu as extensões, a"
        " solução em comprimento de onda, a relação sinal-ruído por ordem e as"
        " velocidades baricêntricas usadas aqui; este pacote parte dos seus"
        " t.fits e nunca volta a derivar nada disso.",
    "help_log":
        "O registo próprio da execução, com as cores que um terminal lhe"
        " daria: verde para o progresso, azul para um número, laranja para algo"
        " ignorado, vermelho para o que para uma execução.",
}
