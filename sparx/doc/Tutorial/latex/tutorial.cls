\ProvidesClass{tutorial}[2016/11/23 v1 wrapper to scrbook for SPHIRE tutorial]

% All other Options go to KOMA-Script
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrbook}}
\PassOptionsToPackage{bibliography=totoc}{biblatex}
\ProcessOptions\relax

\LoadClass[
  fontsize=12pt,
  bibliography=totoc,     % bibliography as unnumbered chapter in toc
  numbers=noenddot,       % no dot after figure/table number
  captions=tableheading,  % correct spacing for table headings
  titlepage=firstiscover, % symmetrical margins on titlepage
  headings=normal         % size of chapter headings slightly smaller
]{scrbook}

% for ogn color
\RequirePackage[hyphens]{url}
% pdf sizes for ogn color
\RequirePackage{luatex85}
% Color
\RequirePackage{xcolor}
% Change chapter/section/...
\RequirePackage[explicit]{titlesec}
% Underline chapter
\RequirePackage{ulem}
% modify the table of contents
\RequirePackage{eso-pic}
% Booleans
\RequirePackage{ifthen}

% Define colors
\xdefinecolor{sphire}{RGB}{103, 162, 195}

\xdefinecolor{chapter_general}{RGB}{23, 54, 93}
\xdefinecolor{chapter_green}{RGB}{91, 129, 138}
\xdefinecolor{chapter_orange}{RGB}{180, 122, 95}
\xdefinecolor{chapter_grey}{RGB}{118, 117, 122}
\xdefinecolor{chapter_brown}{RGB}{133, 118, 106}
\xdefinecolor{chapter_red}{RGB}{123, 77, 65}

\xdefinecolor{section}{RGB}{103, 162, 195}

\xdefinecolor{tip}{RGB}{0, 128, 0}

\xdefinecolor{note}{RGB}{77, 128, 179}

\xdefinecolor{important}{RGB}{255, 0, 0}

\xdefinecolor{advanced}{RGB}{128, 128, 128}
\xdefinecolor{advanced_head}{RGB}{128, 128, 128}

\xdefinecolor{issue}{RGB}{128, 128, 128}
\xdefinecolor{issue_head}{RGB}{255, 0, 0}

\xdefinecolor{red}{RGB}{255, 0, 0}
\xdefinecolor{purple}{RGB}{204, 20, 150}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Pictographs
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newlength{\logowidth} % Variable for the width of the logo in the chapter
\setlength{\logowidth}{2cm} % Set the length of the variable in the chapter

\newcommand{\logoproject}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_project.pdf}}}
\newcommand{\logomovie}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_movie.pdf}}}
\newcommand{\logocter}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_cter.pdf}}}
\newcommand{\logowindow}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_window.pdf}}}
\newcommand{\logoisac}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_isac.pdf}}}
\newcommand{\logoviper}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_viper.pdf}}}
\newcommand{\logomeridien}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_meridien.pdf}}}
\newcommand{\logosort}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_sort3d.pdf}}}
\newcommand{\logolocal}{{\protect\includegraphics[width=\logowidth]{media/sxgui_pictograph_localres.pdf}}}

\newcommand{\chapterlogo}{}
\newcommand{\chapterindent}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Reformat Chapter and section and subsection
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GLOABL

% global settings for all types
\setkomafont{disposition}{\rmfamily\bfseries}

% Underline in SPHIRE colors
\newcommand{\sphireunderline}{\bgroup\markoverwith
{\textcolor{sphire}{\rule[-1ex]{2pt}{1pt}}}\ULon}

%%%% CHAPTER

% Chapter font numberless
\titleformat{name=\chapter,numberless}[display]
{\rmfamily\Huge\color{chapter_general}}
{#1}
{-42pt}
{\sphireunderline{\qquad\qquad\qquad}}

% Chapter font not numberless
\titleformat{name=\chapter}[display]
{\rmfamily\Huge}
{#1}
{-35pt}
{\chapterindent\sphireunderline{\hspace{5cm}}}

% Variable to check if the chapter is used the first time
\newboolean{boolvar}
\setboolean{boolvar}{true}

\newcommand{\mychaptertitle}{}
\newcommand{\mychapter}[3]{
	\ifthenelse{\boolean{boolvar}}{\relax\renewcommand{\chapterlogo}{#3}}{\clearpage\renewcommand{\chapterlogo}{#3}}
	\renewcommand{\mychaptertitle}{\color{#1}#2}
	\begin{minipage}[t]{\logowidth + 1em}
		\vspace{0pt}
		\chapterlogo
	\end{minipage}
	\begin{minipage}[t]{\textwidth - \logowidth - 1em}
		\vspace{0.1\logowidth}
		\chapter{\color{#1}#2}
		\label{#2}
	\end{minipage}
	\setboolean{boolvar}{false}
}

% Space before and after the chapter
\titlespacing*{\chapter}
{0pt}{0cm}{0.5cm}

% Show only chaptername not number in the header
\renewcommand{\sectionmark}[1]{\markright{#1}}

%%%% SECTION

% section font numberless
\titleformat{name=\section,numberless}[block]
{\rmfamily\Large}
{#1}
{0pt}
{}

% section font not numberless
\titleformat{name=\section}[block]
{\rmfamily\Large}
{#1}
{0pt}
{}


\newcommand{\mysection}[1]{{\color{section}\textit{\section{#1}\label{#1}}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Header
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newpagestyle{fancy}{
\setheadrule{0.4pt}{}
\sethead% even
{{\setlength{\logowidth}{1cm}\chapterlogo}}{\textsl{\mychaptertitle}}{}% odd

\setfoot% even
{}{\textsl{\thepage}}{}% odd
}
\pagestyle{fancy}

\newpagestyle{first}{
\setheadrule{0.4pt}{}
\sethead% even
{}{}{\textit{Written by Christos Gatsogiannis with input from the \sphire{} team, format by Markus Stabrin}}% odd

\setfoot% even
{}{\textsl{\thepage}}{}% odd
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Titel page picture
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newlength{\titlepicwidth} % Variable for the width of the logo
\setlength{\titlepicwidth}{25cm} % Set the length of the variable


\newcommand{\TitelPage}{% Name of the new Picture
\put(10,13){% Starting in the lower left corner going up and right
\parbox[b][\titlepicwidth]{1.148\titlepicwidth}{% Box Size
\vfill
\raggedright
\includegraphics[width=\titlepicwidth,height=1.148\titlepicwidth,%height is 1.148 time width
keepaspectratio]{media/titlepage.png}% Logo
\vfill
}}}
