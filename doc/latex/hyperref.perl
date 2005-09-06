# HYPERREF - interface for LaTeX2HTML    -*-perl-*-
# Andreas Hirczy <ahi@itp.tu-graz.ac.at>
#
# install in $LATEX2HTML/styles/

package main;

# Suppress option warning messages:

sub do_hyperref_colorlinks {}
sub do_hyperref_backref {}
sub do_hyperref_pagebackref {}
sub do_hyperref_pdfstartview_FitH {}
sub do_hyperref_breaklinks {}

# low meaning in HTML; ignore silently (by now, TODO)
sub do_cmd_hypertarget {}

#
# Functions
#

#
# \href <url> <text> (converted from htlink.perl)
#
sub do_cmd_href{
    local($_) = @_;
    local($text, $url);
    s/$next_pair_pr_rx/$url  = $2; ''/eo;
    s/$next_pair_pr_rx/$text = $2; ''/eo;
    # and recode the ~ (don't turn it to space)
    $url =~ s/~/&#126;/go;
    join('',"<A HREF=\"$url\">$text</A>",$_);
}


#sub do_cmd_hypersetup{
#}

# Replace `meta_information' in latex2html.config
sub meta_information {
    local($_) = @_;
    # Cannot have nested HTML tags...
    do { s/<[^>]*>//g;
         "<META NAME=\"description\" CONTENT=\"$_\">\n" .
	 "<META NAME=\"resource-type\" CONTENT=\"document\">\n" .
	 "<META NAME=\"distribution\" CONTENT=\"global\">\n$htmetainfo" } if $_;
}


1;                              # This must be the last line
