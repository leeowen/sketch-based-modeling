% ---------------------------------------------------------------------
% EG author guidelines plus sample file for EG publication using LaTeX2e input
% D.Fellner, v2.02, Jan 25, 2017


\title[EG \LaTeX\ Author Guidelines]%
      {Sketch-based character modelling with primitive deformer and detail generator}

% for anonymous conference submission please enter your SUBMISSION ID
% instead of the author's name (and leave the affiliation blank) !!
%\author[D. Fellner \& S. Behnke]
\author[SUBMISSION ID]
{
%\parbox{\textwidth}{\centering D.\,W. Fellner\thanks{Chairman Eurographics Publications Board}$^{1,2}$
%        and S. Behnke$^{2}$ 
%        S. Spencer$^2$\thanks{Chairman Siggraph Publications Board}
%        }
        \\
% For Computer Graphics Forum: Please use the abbreviation of your first name.
%{
%\parbox{\textwidth}{\centering $^1$TU Darmstadt \& Fraunhofer IGD, Germany\\
%         $^2$Graz University of Technology, Institute of Computer Graphics and Knowledge Visualization, Austria
%        $^2$ Another Department to illustrate the use in papers from authors
%             with different affiliations
%       }
%}
}
% ------------------------------------------------------------------------

% if the Editors-in-Chief have given you the data, you may uncomment
% the following five lines and insert it here
%
% \volume{36}   % the volume in which the issue will be published;
% \issue{1}     % the issue number of the publication
% \pStartPage{1}      % set starting page


%-------------------------------------------------------------------------
\begin{document}

% uncomment for using teaser
% \teaser{
%  \includegraphics[width=\linewidth]{eg_new}
%  \centering
%   \caption{New EG Logo}
% \label{fig:teaser}
%}

\maketitle
%-------------------------------------------------------------------------
\begin{abstract}
This paper proposes a new sketch-guided and ODE-driven character modelling technique. Our system consists of two main components: primitive deformer and detail generator. With such a technique, we first draw 2D silhouette contours of a character model. Then, we select proper primitives and align them with the corresponding silhouette contours. After that, we develop a sketch-guided and ODE-driven primitive deformer. It uses ODE-based deformations to deform the primitives to exactly match the generated 2D silhouette contours in one view plane and obtain a base mesh of a character model consisting of deformed primitives. In order to add various 3D details, we develop a local shape generator which uses sketches in different view planes to define a local shape and employs ODE-driven deformations to create a local surface passing through all the sketches. \\
The experimental results demonstrate that our proposed approach can create 3D character models with 3D details from 2D sketches easily, quickly and precisely.
\\
%The ABSTRACT is to be in fully-justified italicized text, 
%   between two horizontal lines,
%in one-column format, 
%   below the author and affiliation information. 
%   Use the word ``Abstract'' as the title, in 9-point Times, boldface type, 
%   left-aligned to the text, initially capitalized. 
%   The abstract is to be in 9-point, single-spaced type.
%   The abstract may be up to 3 inches (7.62 cm) long. \\
%   Leave one blank line after the abstract, 
%   then add the subject categories according to the ACM Classification Index 
%-------------------------------------------------------------------------
%  ACM CCS 1998
%  (see http://www.acm.org/about/class/1998)
% \begin{classification} % according to http:http://www.acm.org/about/class/1998
% \CCScat{Computer Graphics}{I.3.3}{Picture/Image Generation}{Line and curve generation}
% \end{classification}
%-------------------------------------------------------------------------
%  ACM CCS 2012(see http://www.acm.org/about/class/class/2012)
%The tool at \url{http://dl.acm.org/ccs.cfm} can be used to generate
% CCS codes.
%Example:
\begin{CCSXML}
<ccs2012>
<concept>
<concept_id>10010147.10010371.10010396.10010397</concept_id>
<concept_desc>Computing methodologies~Mesh models</concept_desc>
<concept_significance>500</concept_significance>
</concept>
%<concept>
%<concept_id>10010583.10010588.10010559</concept_id>
%<concept_desc>Hardware~Sensors and actuators</concept_desc>
%<concept_significance>300</concept_significance>
%</concept>
</ccs2012>
\end{CCSXML}
\ccsdesc[500]{Computing methodologies~Mesh model}
%\ccsdesc[300]{Hardware~Sensors and actuators}
\printccsdesc   
\end{abstract}  
%-------------------------------------------------------------------------
\section{Introduction}

%Please follow the steps outlined in this document very carefully when
%submitting your manuscript to Eurographics.

%You may as well use the \LaTeX\ source as a template to typeset your own
%paper. In this case we encourage you to also read the \LaTeX\ comments
%embedded in the document.
The main contributions made by our proposed approach are:
\begin{itemize}
  \item[$\bullet$] We develop an efficient sketch-guided and ODE-driven primitive deformer to create a base mesh. It can deform primitives to exactly match the generated silhouette contours. Compared to the existing methods, it automates shape manipulation, avoids tedious manual operations, can deform primitives to match the generated silhouette contours quickly, and is powerful in achieving different shapes of a same primitive.
  \item[$\bullet$] We develop a detail generator to add 3D details to the base mesh. Our proposed sketch-guided and ODE-driven local shape creator can create a new local shape to match user's drawn sketches in different views quickly. The image-based detail generator can automatically generate fine details from 2D images. 
  \item[$\bullet$] Our character modelling system provide editing operations, support users to manipulate the control curves the draw in the creation stage and adding more control curves after the surface has been created. With our developed system, 3D character models with 3D details can be created easily and efficiently.
\end{itemize}
The rest of the paper is organised as follows. The previous related work is briefly reviewed in Section 2. The system overview of our proposed approach is presented in Section 3. Primitive deformer is examined in Section 4, and Detail generator is investigated in Section 5. Finally, the conclusions and future work are discussed in Section 6. 
%-------------------------------------------------------------------------
\section{Related Work}
\begin{description}[wide=0\parindent]
\item[Inflation] 
Over the past two decades, sketch-based-modelling (SBM) has been widely studied in the computer graphic research community. Several research based systems have been proposed to generate organic models. The surface inflation technique extrudes the polygonal mesh from the skeleton outwards do a good job in modelling stuffed toys. One trend is to inflate freeform surfaces to create simple stuffed animals and other rotund objects in a sketch-based modelling fashion, such as \cite{igarashi1999teddy, Karpenko:2006:SFS:1141911.1141928, nealen2007fibermesh}. The Teddy system \cite{igarashi1999teddy} is the pioneer, it takes closed curves as inputs and find their cordial axes as spline, then wrap the splines with the polygonal mesh. Later, FiberMesh\cite{nealen2007fibermesh} enriched the editing operations for the inflating base mesh. FiberMesh also presents two types of the control curves: smooth and sharp. A smooth curve constrains the surface to be smooth across it, while a sharp curve only places positional constraints with C0 continuity. Sharp control curves appears when operations like cutting, extrusion and tunnel take place. Sharp control curves also serves the creation of creases on surface. Based on the study from William\cite{williams1997stochastic}, the SmoothSketch system addressed the problem of T-junction and cusp, which Teddy fails to solve.\\
\item[Primitive] 
Unlike the inflating systems, primitives-based systems deconstruct the modelling task as a process of creating a certain set of geometry primitives and further editing on the primitives. The idea of assembling simple geometry primitives to form 3D models is very common in CSG(constructive solid geometry) modelling, related researches including \cite{shtof2013geosemantic,chen20133}. Shtof et.al \cite{shtof2013geosemantic} introduces a snapping method which helps determining the position and core parameters of several simple geometry primitives. In \cite{chen20133}, the authors provides tools for generating a cylinder from only 3 strokes: the first two strokes define the 2D profile and the last stroke defines the axis along which the profile curve will sweep. Copies of the profile are not only perpendicularly aligned to the axis, but also resized to snap to the input outlines. However their work is only for man-made objects which simple sweeping surface can meet the quality requirements of the shapes. Structured Annotations for 2D-to-3D Modelling \cite{gingold2009structured}, on the other hand, focus on organic modelling. It is a system using two sets of the primitives, one is generalized cylinders, created by the input of a single open sketch stroke representing the spline, and then modified by using simple gestures such as tilt, scale local radius, rotate symmetrical plane, and change cap size; the other primitive is ellipsoid, generated according to the drawn closed ellipse sketch stroke. As the system's name indicates, there are a set pf annotation tools to further editing the surface shape using the annotations such as same-lengths, same angles, alignment and mirror symmetry. \\
\item[Mesh editing] 
A universal disadvantage of above systems is features provided are limited it was not suited to modelling complex and production ready models. This is partly due to their limitations in mesh editing operations, so several sketch-based editing methods have been development to help to improve the performance of sketch-based modelling system. One group of sketch-based editing methods \cite{draper2003gestural, kho2005sketching} treats individual sketch curve as the reference for deformations like bending. Take the \cite{kho2005sketching} for example, user will draw a first curve as a reference curve to both set the ROI (region of interest) and be a controlling 'skeleton', and a new sketch curve indicating the desired deformation of the reference curve. Then mapping the reference curve to the second curve, the ROI will be deformed. These methods are very useful for simple keyframe animation, but for more small-scaled shape's editing task, studies like \cite{lipman2004differential,nealen2007sketch,SBM:SBM07:023-030} are more suitable. These studies take mesh's contours (depend on the viewplan at the time of editing operation is happening) as reference curves, and the newly sketched curves indicating the new positions that the ROI of mesh should be reconstructed to. As they aiming at maintaining the local geometrical details, Laplacian coordinates and correspondingly laplacian surface editing technology are implemented to preserving the features. 
\end{description}

%-------------------------------------------------------------------------
\subsection{Language}

%All manuscripts must be in English.

%-------------------------------------------------------------------------
\subsection{Margins and page numbering}

All printed material, including text, illustrations, and charts,
must be kept within a print area %6.31 inches (16.03 cm) 
7 inches (17.7 cm) wide by
%9.10 inches (23.13 cm) 
9.44 inches (24 cm) high. Do not write or print anything
outside the print area. Number your pages on odd sites right
above, on even sites left above, no page number on the first site.
Do not use page numbering within the final version of your paper.


%------------------------------------------------------------------------
\subsection{Formatting your paper}

% All text with the exception of the abstract must be in a two-column format.
% The total allowable width of the text area -- including header and footer
% lines -- is 161\,mm (6.34 inch) wide by 231\,mm (9.10 inch) high.
% 
% Columns are to be 76\,mm (3.0 inch) wide, with a 8\,mm (0.315 inch) space 
% between them.

All text with the exception of the abstract must be in a two-column format.
The total allowable width of the text area -- including header and footer
lines -- is 177\,mm (7 inch) wide by 245\,mm (9.64 inch) high.

Columns are to be 84\,mm (3.3 inch) wide, with a 8\,mm (0.315 inch) space 
between them.

The space between the header line and the first line of the text body and
between the last line of the text body and the footer line is 5\,mm
(0.196 inch) each.

%-------------------------------------------------------------------------
\subsection{Type-style and fonts}

Wherever Times is specified, Times Roman may also be used. If
neither is available on your word processor, please use the font
closest in appearance to Times that you have access to. Only
Type-1 fonts will be accepted.

MAIN TITLE. The title should be in Times 17-point, boldface type and
centered. Capitalize the first letter of nouns, pronouns, verbs, adjectives,
and adverbs; do not capitalize articles, coordinate conjunctions, or
prepositions (unless the title begins with such a word). Leave two blank
lines after the title.

AUTHOR NAME(s) and AFFILIATION(s) are to be centered beneath the title and
printed in Times 9-point, non-boldface type. This information is to be
followed by two blank lines.

The ABSTRACT ist to be in a one-column format. The MAIN TEXT is to be in a
two-column format.

MAIN TEXT. Type main text in 9-point Times, single-spaced. Do \emph{not} use
double-spacing. All paragraphs should be indented 1 em (the length of the
dash in the actual font). Make sure your text is fully justified -- that is,
flush left and flush right. Please do not place any additional blank lines
between paragraphs. Figure and table captions should be 9-point Times
boldface type as in Figure~\ref{fig:firstExample}.

\noindent Long captions should be set as in Figure~\ref{fig:ex1} or
Figure~\ref{fig:ex3}.

\begin{figure}[htb]
   % an empty figure just consisting of the caption lines
   \caption{\label{fig:ex1}
     'Empty' figure only to serve as an example of long caption requiring 
     more than one line. It is not typed centered but aligned on both sides.}
\end{figure}

\noindent
Figures which need the full textwidth can be typeset as Figure~\ref{fig:ex3}.

\noindent Callouts should be 9-point Times, non-boldface type. Initially
capitalize only the first word of section titles and first-, second-, and
third-order headings.

FIRST-ORDER HEADINGS. (For example, \textbf{1. Introduction}) should be Times
9-point boldface, initially capitalized, flush left, with one blank line
before, and one blank line after.

SECOND-ORDER HEADINGS. (For example, \textbf{2.1. Language}) should be Times
9-point boldface, initially capitalized, flush left, with one blank line
before, and one after. If you require a third-order heading (we discourage
it), use 9-point Times, boldface, initially capitalized, flush left, preceded
by one blank line, followed by a period and your text on the same line.

The headline \emph{(authors / title)} must be shortened if it uses the full 
two column width of the main text.
There must be enough space for the page numbers. Please use ``et al.'' if 
there are more than three authors and specify a shortened version for your title.
%-------------------------------------------------------------------------
\subsection{Footnotes}

Please do \emph{not} use footnotes at all!


%-------------------------------------------------------------------------
\subsection{References}

List all bibliographical references in 9-point Times, single-spaced, at the
end of your paper in alphabetical order. When referenced in the text, enclose
the citation index in square brackets, for example~\cite{Lous90}. Where
appropriate, include the name(s) of editors of referenced books.

For your references please use the following algorithm:
\begin{itemize} 
\item \textbf{one} author: first 3 chars plus year -- 
      e.g.\ \cite{Lous90}
\item \textbf{two}, \textbf{three} or \textbf{four} authors: first char
      of each family name plus year --  e.g.\ \cite{Fellner-Helmberg93} 
      or \cite{Kobbelt97-USHDR} or \cite{Lafortune97-NARF}
\item \textbf{more than 4} authors: first char of family name from 
      first 3 authors followed by a '*' followed by the year -- 
      e.g.\ \cite{Buhmann:1998:DCQ} or \cite{FolDamFeiHug.etal93} 
\end{itemize}

For BibTeX users a style file \ \texttt{eg-alpha.bst} \ is available which
uses the above algorithm.

%-------------------------------------------------------------------------
\subsection{Illustrations, graphs, and photographs}

All graphics should be centered.

%%%
%%% Figure 1
%%%
\begin{figure}[htb]
  \centering
  % the following command controls the width of the embedded PS file
  % (relative to the width of the current column)
  \includegraphics[width=.8\linewidth]{sampleFig}
  % replacing the above command with the one below will explicitly set
  % the bounding box of the PS figure to the rectangle (xl,yl),(xh,yh).
  % It will also prevent LaTeX from reading the PS file to determine
  % the bounding box (i.e., it will speed up the compilation process)
  % \includegraphics[width=.95\linewidth, bb=39 696 126 756]{sampleFig}
  %
  \parbox[t]{.9\columnwidth}{\relax
           For all figures please keep in mind that you \textbf{must not}
           use images with transparent background! 
           }
  %
  \caption{\label{fig:firstExample}
           Here is a sample figure.}
\end{figure}

If your paper includes images, it is very important that they are of
sufficient resolution to be faithfully reproduced.

To determine the optimum size (width and height) of an image, measure
the image's size as it appears in your document (in millimeters), and
then multiply those two values by 12. The resulting values are the
optimum $x$ and $y$ resolution, in pixels, of the image. Image quality
will suffer if these guidelines are not followed.

Example 1: 
%
An image measures 50\,mm by 75\,mm when placed in a document. This
image should have a resolution of no less than 600 pixels by 900
pixels in order to be reproduced faithfully.

Example 2: 
%
Capturing a screenshot of your entire $1024 \times 768$ pixel display
monitor may be useful in illustrating a concept from your research. In
order to be reproduced faithfully, that $1024 \times 768$ image should
be no larger than 85 mm by 64 mm (approximately) when placed in your
document.


%-------------------------------------------------------------------------
\subsection{Color}

\textbf{Please observe:} as of 2003 publications in the proceedings of the
Eurographics Conference can use color images throughout the paper. No
separate color tables are necessary.

However, workshop proceedings might have different agreements! 
Figure~\ref{fig:ex3} is an example for creating color plates.

%------------------------------------------------------------------------
\subsection{Embedding of Hyperlinks / Typesetting of URLs}

Due to the use of the package \texttt{hyperref} the original behavior
of the command $\backslash$\texttt{url} from the package \texttt{url}
is not available. To circumvent this problem we either recommend to
use the command $\backslash$\texttt{httpAddr} from the 
included package \texttt{egweblnk} (see below) or to replace the
command $\backslash$\texttt{url} by the command $\backslash$\texttt{webLink} 
-- e.g. in cases where $\backslash$\texttt{url} has been used
widely in BibTeX-References. In the latter case we suggest to run
BibTeX as usual and then replace all occurences of $\backslash$\texttt{url}  by
$\backslash$\texttt{webLink}

\noindent
The provided commands for hyperlinks, in a nutshell, are:

\begin{description} \itemsep 1ex
\item [\webLinkFont $\backslash$httpAddr \{URL without leading 'http:'\}]
      \mbox{}\\
      e.g. \  \httpAddr{//diglib.eg.org/EG/DL/WS}

\item[\webLinkFont $\backslash$ftpAddr \{URL without leading 'ftp:'\}]
      \mbox{}\\
      e.g. \  \ftpAddr{//www.eg.org/EG/DL/ftpupload}

\item[\webLinkFont $\backslash$URL \{url\}]
      \mbox{}\\
      e.g. \  \URL{http://www.eg.org/EG/DL/WS}

\item[\webLinkFont $\backslash$MailTo \{Email addr\}]
      \mbox{}\\
      e.g. \  \MailTo{publishing@eg.org}

\item[\webLinkFont $\backslash$MailToNA \{emailName\}\{@emailSiteAddress\}]
      \mbox{}\\
      e.g. \  \MailToNA{publishing}{@eg.org}

\item[\webLinkFont $\backslash$webLink\{URL without hyperlink creation\}]
      \mbox{}\\
      e.g. \  \webLink{http://www.eg.org/some_arbitrary_long/but_useless/URL}

\end{description}


%------------------------------------------------------------------------
\subsection{Comparison with XXX}

Our approach is compared against XXXX, which is the only method for XXX whose code is freely available. 
Its main disadvantage is that it is not suited for production ready models that can be modelled with professional modelling packages such Maya or ZBrush.\\
%------------------------------------------------------------------------
\subsection{Future work}

For now in the modelling stage of our system, all our control curves are creating smooth surfaces with C2 continuity, and the result template mesh looks pleasing as it is. In our future work, we will investigate the potentials of creases in more detailed character sculpting.

Your final paper should be delivered as a PDF document with all typefaces
embedded. \LaTeX{} users should use \texttt{dvips} and \texttt{ps2pdf} to
create this PDF document. Adobe Acrobat Distiller may be used in place of
\texttt{ps2pdf}.

Adobe PDFWriter is \emph{not} acceptable for use. Documents created with
PDFWriter will be returned to the author for revision. \texttt{pdftex} and
\texttt{pdflatex} (and its variants) can be used only if the author can
make certain that all typefaces are embedded and images are not downsampled
or subsampled during the PDF creation process.

Users with no access to these PDF creation tools should make available a
PostScript file and we will make a PDF document from it.


The PDF file \emph{must not} be change protected.


%------------------------------------------------------------------------
\subsubsection*{Configuration Notes: dvips / ps2pdf / etc.}

\noindent
\texttt{dvips} should be invoked with the \texttt{-Ppdf} and \texttt{-G0}
flags in order to use Type 1 PostScript typefaces:

\begin{verbatim}
    dvips -t a4 -Ppdf -G0 -o my.ps my.dvi
\end{verbatim}


\noindent
If you are using version 7.x of GhostScript, please use the following method of invoking \texttt{ps2pdf}, in
order to embed all typefaces and ensure that images are not downsampled or subsampled in the PDF
creation process:

\begin{verbatim}
  ps2pdf -dMaxSubsetPct=100 \
         -dCompatibilityLevel=1.3 \
         -dSubsetFonts=true \
         -dEmbedAllFonts=true \
         -dAutoFilterColorImages=false \
         -dAutoFilterGrayImages=false \
         -dColorImageFilter=/FlateEncode \
         -dGrayImageFilter=/FlateEncode \
         -dMonoImageFilter=/FlateEncode \
         mypaper.ps mypaper.pdf
\end{verbatim}


If you are using version 8.x of GhostScript, please use this method in place of the example above:
\begin{verbatim}
  ps2pdf -dPDFSETTINGS=/prepress \
         -dCompatibilityLevel=1.3 \
         -dAutoFilterColorImages=false \
         -dAutoFilterGrayImages=false \
         -dColorImageFilter=/FlateEncode \
         -dGrayImageFilter=/FlateEncode \
         -dMonoImageFilter=/FlateEncode \
         -dDownsampleColorImages=false \
         -dDownsampleGrayImages=false \
         mypaper.ps mypaper.pdf
\end{verbatim}

%------------------------------------------------------------------------
\subsubsection*{Configuration Notes: pdftex / pdflatex / etc.}

\noindent
Configuration of these tools to embed all typefaces can be accomplished by editing the \texttt{updmap.cfg} file
to enable inclusion of the standard (or base) 14 typefaces.

Linux users can run the \texttt{updmap} script to do this:
\begin{verbatim}
updmap --setoption pdftexDownloadBase14 true
\end{verbatim}

Windows users should edit the \texttt{updmap.cfg} files found in their TeX installation directories (one or both
of the following may be present):
\begin{verbatim}
  INSTALLDIR\texmf\web2c\updmap.cfg
  INSTALLDIR\localtexmf\miktex\config\updmap.cfg
\end{verbatim}

Ensure the value for \texttt{pdftexDownloadBase14} is "true," and then follow the instructions found here:
\httpAddr{//docs.miktex.org/manual/} to update your MikTeX installation.

%------------------------------------------------------------------------
\subsubsection*{Configuration Notes: Acrobat Distiller}

We recommend to use a Distiller job options file that embeds
all typefaces and does not downsample or subsample images when creating the PDF document.
%------------------------------------------------------------------------
\subsection{Exclusive License Form}

You must include your signed Eurographics Exclusive License Form 
when you submit your finished paper. We MUST have this form before
your paper can be published in the proceedings.

%-------------------------------------------------------------------------
\subsection{Conclusions}
%Please direct any questions to the production editor in charge of
%these proceedings.
Our proposed system on the other hand aims at creating models representing various types of characters with comparatively more detail with the same level of effort and in the same time. Based on the , while also aims at enhancing it by: 1. Introduce superquadric primitives to better convey the organic feeling; 2. precisely fitting the simple geometrical primitives' contour into the input sketch strokes while adjusting the overall shape aesthetically and physical plausibly with the help of a new ODE deforming technology; 3. generate patches which supports multiple numbers of control curves, in order to achieve better quality models with better accuracy and details and 4. Moreover this system support the users to modify the shape by editing the sketch strokes. 
%-------------------------------------------------------------------------

%\bibliographystyle{eg-alpha}
\bibliographystyle{eg-alpha-doi}

\bibliography{egbibsample}

%-------------------------------------------------------------------------
\newpage


\begin{figure*}[tbp]
  \centering
  \mbox{} \hfill
  % the following command controls the width of the embedded PS file
  % (relative to the width of the current column)
  \includegraphics[width=.3\linewidth]{sampleFig}
  % replacing the above command with the one below will explicitly set
  % the bounding box of the PS figure to the rectangle (xl,yl),(xh,yh).
  % It will also prevent LaTeX from reading the PS file to determine
  % the bounding box (i.e., it will speed up the compilation process)
  % \includegraphics[width=.3\linewidth, bb=39 696 126 756]{sampleFig}
  \hfill
  \includegraphics[width=.3\linewidth]{sampleFig}
  \hfill \mbox{}
  \caption{\label{fig:ex3}%
           For publications with color tables (i.e., publications not offering
           color throughout the paper) please \textbf{observe}: 
           for the printed version -- and ONLY for the printed
           version -- color figures have to be placed in the last page.
           \newline
           For the electronic version, which will be converted to PDF before
           making it available electronically, the color images should be
           embedded within the document. Optionally, other multimedia
           material may be attached to the electronic version. }
\end{figure*}

\end{document}
