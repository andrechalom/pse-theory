R --vanilla <<EOF
	library(utils)
	Sweave("$1.Rnw", encoding="utf8")
EOF
