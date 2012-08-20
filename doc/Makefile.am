MANNAME = manual
MANSRC = $(MANNAME).tex
MANPDF = $(MANNAME).pdf

if HAVE_PDFLATEX

doc_DATA = $(MANPDF)

$(MANPDF): $(srcdir)/$(MANSRC)
	    $(PDFLATEX) $<
	    $(PDFLATEX) $<

CLEANFILES = \
	$(MANNAME).aux \
	$(MANNAME).log \
	$(MANNAME).toc \
	$(MANPDF)

endif

EXTRA_DIST = \
	$(MANSRC)

MAINTAINERCLEANFILES = \
	Makefile.in
