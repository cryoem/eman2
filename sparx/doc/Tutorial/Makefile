TEXOPTIONS = -lualatex \
			 --output-directory=../build \
			 --interaction=nonstopmode \
			 --halt-on-error

TEXGLOBAL = latex/tutorial.tex \
			latex/tutorial.cls \
			latex/Tex_global/header_1.tex \
			latex/Tex_global/header_2.tex \
			latex/Tex_global/header_font.tex \
			latex/Tex_global/lit.bib \
			latex/Tex_global/acronyms.tex

TEXCONTENT = build/tutorial_parsed.tex

CURRENT_TIME = $(shell date '+%Y_%m_%d_%H_%M_%S')

all: build/tutorial.pdf

build/tutorial.pdf: $(TEXGLOBAL) $(TEXCONTENT) | build
	cd latex; latexmk $(TEXOPTIONS) tutorial.tex >> ../log.txt 2>> ../log.txt

$(TEXCONTENT): tutorial.docx | build versions
	python scripts/docx_to_gfm.py tutorial.docx versions/tutorial_$(CURRENT_TIME).gfm --tags_file scripts/gfm_to_latex_tags.txt >> log.txt 2>> log.txt
	python scripts/gfm_to_latex.py versions/tutorial_$(CURRENT_TIME).gfm $(TEXCONTENT) --acronyms_file latex/Tex_global/acronyms.tex --tags_file scripts/gfm_to_latex_tags.txt >> log.txt 2>> log.txt

latest_to_docx:
	$(eval time_stemp=$(shell ls -r versions | awk 'NR==1{print$1}' | sed 's/tutorial_\(.*\).gfm/\1/g'))
	python scripts/gfm_to_docx.py versions/tutorial_${time_stemp}.gfm tutorial_${time_stemp}.docx

check_docx:
	python scripts/docx_to_gfm.py tutorial.docx tutorial_check.gfm --tags_file scripts/gfm_to_latex_tags.txt

readme:
	pandoc README.md -o README.pdf

versions:
	mkdir -p versions

build:
	mkdir -p build

clean:
	rm -rf build
	rm -rf log.txt
