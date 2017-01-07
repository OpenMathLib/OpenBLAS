#!/bin/bash

while read OBJ; do
	
	if echo "$OBJ"|grep "_$" >/dev/null
	then

		O1=$(echo "$OBJ"|sed -e 's/_$//' )

		if grep -w "$O1" exports/gensymbol >/dev/null
        	then
			true
		else
			echo "$O1"
		fi
	fi

	if echo "$OBJ"|grep "^cblas" >/dev/null
	then

		if grep -w "$OBJ" exports/gensymbol >/dev/null
        	then
			true
		else
			echo "$OBJ"
		fi
	fi

	if echo "$OBJ"|grep "^LAPACKE" >/dev/null
	then

		if grep -w "$OBJ" exports/gensymbol >/dev/null
        	then
			true
		else
			echo "$OBJ"
		fi
	fi

	if echo "$OBJ"|grep "^lapack" >/dev/null
	then

		if grep -w "$OBJ" exports/gensymbol >/dev/null
        	then
			true
		else
			echo "$OBJ"
		fi
	fi




done

