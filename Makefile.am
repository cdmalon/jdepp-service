SUBDIRS = src tools
EXTRA_DIST = GPL LGPL tools

model_topdir = model
corpus = @CORPUS@
modeldir = ${model_topdir}/${corpus}
CLEANFILES = ${modeldir}/*
JDEPP = src/jdepp$(EXEEXT)
KNBC_URL = http://nlp.ist.i.kyoto-u.ac.jp/kuntt/KNBC_v1.0_090925.tar.bz2

CHUNKER_LEARNER_OPT=-t 1 -d 2 -c 0.001 -i 40
CHUNKER_CLASSIFIER_OPT=-s 0.03 -i 14
PARSER_LEARNER_OPT=-t 1 -d 2 -c 0.00005 -i 40
PARSER_CLASSIFIER_OPT=-s 0.015 -i 8

# fastest configuration without conjunctive features
# CHUNKER_LEARNER_OPT="-t 0 -c 0.05 -i 40"
# CHUNKER_CLASSIFIER_OPT=""
# PARSER_LEARNER_OPT="-t 0 -c 0.001 -i 40"
# PARSER_CLASSIFIER_OPT=""

install-data-hook:
	@echo "Installing ${corpus} model.."
	if test -e ${modeldir}/dic.euc; then \
		$(INSTALL) -d $(DESTDIR)$(pkglibdir)/${modeldir}; \
		$(INSTALL) -p ${modeldir}/* $(DESTDIR)$(pkglibdir)/${modeldir}; \
	echo 'done.'; \
	fi
	@if ! test -e ${modeldir}/dic.euc; then \
		echo 'WARNING: No model has been built; try make model.'; \
	fi

uninstall-local:
	if test -e $(DESTDIR)$(pkglibdir)/${modeldir}; then \
		$(RM) $(DESTDIR)$(pkglibdir)/${modeldir}/*; \
	fi

knbc:
	$(MKDIR_P) ${model_topdir}/knbc
	$(RM) ${model_topdir}/knbc/*
	@if ! test -e "$(CORPUS_DIR)"; then \
		echo "Kyoto-University and NTT Blog (KNB) Corpus not found: $(CORPUS_DIR)"; \
		echo "Downloading KNB Corpus.."; \
		if test $(CURL) != "no"; then \
			$(CURL) -O $(KNBC_URL); \
		elif test $(WGET) != "no"; then \
			$(WGET) $(KNBC_URL); \
		else \
			echo "no downloaders (curl, wget); I can't do anything more"; \
			exit 1; \
		fi; \
		echo "done."; \
		echo "Extracting tarball.."; \
		tar jxf KNBC_v1.0_090925.tar.bz2; \
		echo "done."; \
	fi
	find $(CORPUS_DIR)/corpus1 -type f -name "KN*" | LC_ALL=C sort | $(AWK) -F'N|_' '$$(NF-2) % 6 <  4 { print }' | xargs cat | @KNBC_CONV@ > ${model_topdir}/knbc/train.JDP.euc
	find $(CORPUS_DIR)/corpus1 -type f -name "KN*" | LC_ALL=C sort | $(AWK) -F'N|_' '$$(NF-2) % 6 == 4 { print }' | xargs cat | @KNBC_CONV@ > ${model_topdir}/knbc/dev.JDP.euc
	find $(CORPUS_DIR)/corpus1 -type f -name "KN*" | LC_ALL=C sort | $(AWK) -F'N|_' '$$(NF-2) % 6 == 5 { print }' | xargs cat | @KNBC_CONV@ > ${model_topdir}/knbc/test.JDP.euc
	@echo "Converting corpus into JDP format.."

kyoto:
	$(MKDIR_P) ${model_topdir}/kyoto
	$(RM) ${model_topdir}/kyoto/*
	@test -e "$(CORPUS_DIR)" || (echo 'Kyoto University Text Corpus not found (CORPUS_DIR=$(CORPUS_DIR))' && exit 1)
	@echo "Converting corpus into JDP format.."
	cat $(CORPUS_DIR)/dat/syn/95{01{01,03,04,05,06,07,08,09,10,11},{01,02,03,04,05,06,07,08}ED}.KNP | @KYOTO_CONV@ > ${model_topdir}/kyoto/train.JDP.euc
	cat $(CORPUS_DIR)/dat/syn/95{01{12,13},09ED}.KNP | @KYOTO_CONV@ > ${model_topdir}/kyoto/dev.JDP.euc
	cat $(CORPUS_DIR)/dat/syn/95{01{14,15,16,17},{10,11,12}ED}.KNP | @KYOTO_CONV@ > ${model_topdir}/kyoto/test.JDP.euc
	@echo "done.\n"

${modeldir}/dic.euc: model

model: all ${corpus}
	@echo "translate corpus into UTF-8.."
	for i in train dev test; do \
	  $(ICONV) -f EUC-JP -t UTF-8 ${modeldir}/$$i.JDP.euc > ${modeldir}/$$i.JDP.utf8; \
	done
	@echo "done."
	@echo "Training and testing BUNSETSU CHUNKER.."
	${JDEPP} -t 0 -I 1 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${CHUNKER_LEARNER_OPT}
	${JDEPP} -t 3 -I 1 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${CHUNKER_LEARNER_OPT} -- ${CHUNKER_CLASSIFIER_OPT} -t 0
	${JDEPP} -t 1 -I 1 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${CHUNKER_LEARNER_OPT} -- ${CHUNKER_CLASSIFIER_OPT} < ${modeldir}/test.JDP.utf8
	@echo "done."
	@echo "Training and testing DEPENDENCY PARSER.."
	${JDEPP} -t 0 -I 2 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${PARSER_LEARNER_OPT}
	${JDEPP} -t 3 -I 2 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${PARSER_LEARNER_OPT} -- -- ${PARSER_CLASSIFIER_OPT} -t 0
	${JDEPP} -t 1 -I 2 -c ${modeldir}/train.JDP.utf8 -m ${modeldir} -- ${PARSER_LEARNER_OPT} -- -- ${PARSER_CLASSIFIER_OPT} < ${modeldir}/test.JDP.utf8
	@echo "done."
	@echo "Generating feature strings in EUC-JP"
	${JDEPP} -t 1 -e 1 -I 1 -c ${modeldir}/train.JDP.euc -m ${modeldir} < ${modeldir}/test.JDP.euc 2> /dev/null
	@echo "done."

.PHONY: model knbc kyoto
