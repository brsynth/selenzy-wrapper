function customConfirm( message, title ) {
    console.log('here');
    if (!title ) 
	title = 'Alert';
    if ( !message ) 
	message = 'System error';
    $( '<div></div>' ).html( message ).dialog( {
	title: title,
	resizable: false,
	modal: true,
	buttons: {
	    'Confirm': function() {
		deleteRows();
		$( this ).dialog( 'close' )
	    },
	    'Cancel': function() {
		$( this ).dialog( 'close' )
	    }
	   }
	});
}


function formatTable(csvdata) {
    $('.dataframe').addClass('Selenzy');
    $('.Selenzy th').each( function(index, value) {
	$( this ).addClass('active');
	});
    addSorted();
    addLinks();
    addSelection();
    if (totalRows() == 0) {
	// If no row is left, disable everything
	$( '.Custom' ).addClass('disabledbutton');
	$( '.Navigate' ).addClass('disabledbutton');
	$( '.Filter' ).addClass('disabledbutton');
    } else {
	$( '.Custom' ).removeClass('disabledbutton');
	$( '.Navigate' ).removeClass('disabledbutton');
	$( '.Filter' ).removeClass('disabledbutton');
    }
}

function formatHeader(filt) {
    for (i=0; i < filt.length; i++) {
	xval = filt[i];
	if (xval > 0) {
	    $( $('thead th')[xval]).attr("sort", "up").attr("title", "Sort by column");
	   }
    }
}


function sortTable(valueSelected) {
    // post and retrieve data from server to sort table
    $.ajax({	
	data : {
	    filter : JSON.stringify(valueSelected),
	    csv : JSON.stringify(csvlink),
	    session: JSON.stringify(sessionid),
	    event: JSON.stringify(event.timeStamp)
	},
	type : 'POST',
	url : '/sorter',
	success : function(serverdata) {
	    data = JSON.parse(serverdata);
	    $('.Selenzy').replaceWith(data.data.csv);
	    formatTable();
	    formatHeader(data.data.filter);
	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}

function addRows() {
    // Add rows from FASTA file
     $( '#add' ).change( function() {
	 if ( $( '#fasta' ).val() != '' ) {
	     $( '#add' ).submit();
	 }
     });

    $( '#add' ).submit(
	function( event ) {
	    var mes = $( '.labfile' ).text();
	    $( '.labfile' ).text('[Submitting sequences..]');
	    $.ajax({
		data: new FormData( this ),
		processData: false,
		contentType: false,		
		type : 'POST',
		url : '/adder',
		success : function(serverdata) {
		    data = JSON.parse(serverdata);
		    $('.Selenzy').replaceWith(data.data.csv);
		    formatTable();
		    $( '.labfile' ).text(mes);
		    // Avoid propagation of the click event
		    $( '#fasta' ).val('');
		    $( '.Remove' ).addClass('disabledbutton');
		},
		error : function() {
		    console.log('Sorry, no luck');
		}
	    });
	});
}



function deleteRows() {
    // Delete selected rows
    $.ajax({	
	data : {
	    filter : JSON.stringify(selectedRows()),
	    csv : JSON.stringify(csvlink),
	    session: JSON.stringify(sessionid),
	    event: JSON.stringify(event.timeStamp)
	},
	type : 'POST',
	url : '/remover',
	success : function(serverdata) {
	    data = JSON.parse(serverdata);
	    $('.Selenzy').replaceWith(data.data.csv);
	    // Avoid propagation of the click event
	    formatTable();
	    $( '.Remove' ).addClass('disabledbutton');

	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}

function updateScore() {
    // Update score
    $.ajax({	
	data : {
	    score : JSON.stringify(readScore().score),
	    csv : JSON.stringify(csvlink),
	    session: JSON.stringify(sessionid),
	    event: JSON.stringify(event.timeStamp)
	},
	type : 'POST',
	url : '/scorer',
	success : function(serverdata) {
	    data = JSON.parse(serverdata);
	    $('.Selenzy').replaceWith(data.data.csv);
	    // Avoid propagation of the click event
	    formatTable();

	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}


function selectedRows() {
    // Returns an array with the indices of selected rows
    var index = [];
    $.each( $( 'tbody th.selected'), function( i, value) {
	index[i] = $( value ).text();
    });
    return(index);

}

function totalRows() {
    // Count number of available rows
    return( $( 'tbody th' ).length );
}


function addColID() {
       $('thead th').each(
	   function( index, value) {
	       if (index > 0) {
		   var title = $( this ).html();
		   $( this ).html( String.fromCharCode( 65 + index - 1) +". "+title );
		  }
	   }
	   );
}


function addSorted() {
       $('thead th').each(
	   function( index, value) {
	       $( this ).attr( "index", index.toString());
	       $( this ).attr( "sort", "down");
	       $( this ).attr( "title", "Sort by column");
	       $( this ).click(
		   function() {
		       if (this.getAttribute("sort") == "up") {
			   $( this ).attr("sort", "down");
			   sortTable([-parseInt(this.getAttribute("index"))]);
		      } else {
			   $( this ).attr("sort", "up");
			   sortTable([parseInt(this.getAttribute("index"))]);
		      }
		   }
		   );
	   }
	   );
}

function addNavigation() {
    var csvTag = $( '<a>' ).attr('href', csvlink)
	.addClass('item').append( "[Download CSV]");
    var fastaTag = $( '<a>' ).attr('href', fastalink).attr('target', '_blank')
	.addClass('item').append("[Download FASTA]");
    var msaTag = $( '<a>' ).attr('href', msafastalink).attr('target', '_blank')
	.addClass('item msa').append("[Download MSA]");
    var msaViewTag = $( '<a>' ).attr('href', msaviewlink).attr('target', '_blank')
	.addClass('item msa').append("[View MSA]");
    if (flagFasta == "False") {
	fastaTag.addClass('disabledbutton');
    }
    if (flagMSA == "False") {
	msaTag.addClass('disabledbutton');
	msaViewTag.addClass('disabledbutton');
    }

    $( '.Navigate' ).append(csvTag).append(fastaTag).append(msaTag).append(msaViewTag);

    addRows();
    $( '.Remove' ).addClass('disabledbutton');
    $( '.Remove').click( function(event) {
	if (selectedRows().length > 0) {
	    customConfirm('Do you want to remove permanently selected rows?');

	  }
    });

}


function addLinks() {
    var seqidix = 1;
    var rxnidix = 5;
    var rows = document.getElementsByTagName('tr');
    for (var i = 1; i < rows.length; i++)	{

	var seqID = rows[i].getElementsByTagName('td')[seqidix];
	var rxnID = rows[i].getElementsByTagName('td')[rxnidix];

	var seqlink = "//www.uniprot.org/uniprot/" + rows[i].getElementsByTagName('td')[seqidix].innerHTML;

	var rxnlink = "http://www.metanetx.org/cgi-bin/mnxweb/equa_info?equa=" + rows[i].getElementsByTagName('td')[rxnidix].innerHTML;

	var aTag = document.createElement('a');
	var bTag = document.createElement('a');
	

	aTag.setAttribute('href', seqlink);
	aTag.setAttribute('target', '_blank');
	bTag.setAttribute('href', rxnlink);
	bTag.setAttribute('target', '_blank');

	aTag.innerHTML = rows[i].getElementsByTagName('td')[seqidix].innerHTML;
	bTag.innerHTML = rows[i].getElementsByTagName('td')[rxnidix].innerHTML;

	seqID.replaceChild(aTag, seqID.childNodes[0]);
	rxnID.replaceChild(bTag, rxnID.childNodes[0]);

    }
}


function addSelection() {
    $('.Selenzy tr th').click(
	function(event) {
	    if ($( this ).hasClass('selected')) {
		$( this ).removeClass('selected');
		$( this ).parent().removeClass('selected');

		if (event.shiftKey) {
		    document.getSelection().removeAllRanges();
		    $( this ).addClass('currentunselection');
		    var last = -1;
		    var current = -1;
		    var rows = [];
		    $.each( $( 'tbody th'), function( i, value) {
			rows[i] = $( value );
			if ( $( value ).hasClass('currentunselection') ) {
			    current = i;
			    } 
			if ( $( value).parent().hasClass('lastunselected') ) {
			    last = i;
			    }
		    });
		    if ( (last >= 0) & (current >= 0) ) {
			for ( i= Math.min(last, current); i != Math.max(last, current); i++ ) {
			    rows[i].removeClass('selected');
			    rows[i].parent().removeClass('selected');
			}
		    }

		    $( this ).removeClass('currentunselection');
		}



		if ( selectedRows().length == 0 ) {
		    $( '.Remove' ).addClass('disabledbutton');
		}
		$.each( $( 'tbody th'), function( i, value) {
		    $( value ).parent().removeClass('lastselected');
		    $( value ).parent().removeClass('lastunselected');
		});
		$( this ).parent().addClass('lastunselected');

	    } else {
		$( this ).addClass('selected');
		$( this ).parent().addClass('selected');			
		$( '.Remove' ).removeClass('disabledbutton');

		if (event.shiftKey) {
		    document.getSelection().removeAllRanges();
		    $( this ).addClass('currentselection');
		    var last = -1;
		    var current = -1;
		    var rows = [];
		    $.each( $( 'tbody th'), function( i, value) {
			rows[i] = $( value );
			if ( $( value ).hasClass('currentselection') ) {
			    current = i;
			    } 
			if ( $( value).parent().hasClass('lastselected') ) {
			    last = i;
			    }
		    });
		    if ( (last >= 0) & (current >= 0) ) {
			for ( i= Math.min(last, current); i != Math.max(last, current); i++ ) {
			    rows[i].addClass('selected');
			    rows[i].parent().addClass('selected');
			}
		    }

		    $( this ).removeClass('currentselection');
		}

		$.each( $( 'tbody th'), function( i, value) {
		    $( value ).parent().removeClass('lastselected');
		    $( value ).parent().removeClass('lastunselected');
		});
		$( this ).parent().addClass('lastselected');

	    }
	}
    );
    
    $('.Selenzy tbody tr th').attr('title', 'Select row');
    lastSelection = 0;
}

function readScore() {
    score = [];
    equation = "<b>Score</b> = ";
    head = true;
    $('.score li').each( function(index, value) {
		if ($( this ).children('.sccheck')[0].checked == true) {
		    val = parseFloat($( this ).children('.scval')[0].value);
		    coln = $( this ).children('.scval').attr('colname');
		    if (val != 0) {
			if (val > 0 ) {
			    equation += '+' ;
			}
			equation +=  val;
			equation +=  ' <b>' + coln + '</b> ';
			}
		    score.push( [coln, val] );
		}
	    });
   return {equation: equation, score: score} ;
}

function initScore() {
    $( '.Filter').click( function(event) {
	$('.score').removeAttr('visibility');
	$('.score').dialog({
		title: 'Sequence score',
	    width: 600,
		height: 400,
	});
    });

	    
    $('.score li').each( function(index, value) {
	if ($( this ).children('.sccheck')[0].checked == false) {
	    $( this ).children('.scval').attr('disabled', 'disabled');
	}
    });
    
    $('.score :checkbox').change( function() {
	if ($(this)[0].checked == true) {
	    $(this).parent().children('.scval').removeAttr('disabled', 'disabled');
	} else {
	    $(this).parent().children('.scval').attr('disabled', 'disabled');
	}
	$('.equation').html( readScore().equation );
    });

    $(".score input[type='number']").bind('click keyup', function() {
	$('.equation').html( readScore().equation );
    });


    $('.equation').html( readScore().equation );

    $('.sccheck').click( function(event) {
	var scolen = readScore().score.length;
	var checkbox = $(this);
	if ((scolen==0) & (!checkbox.is(":checked"))) {
	    checkbox.attr("checked", true)
	    event.preventDefault();
	    return false;
	}
    });

    $( '.updatescore').click( function(event) {
	updateScore();
    });

}

// Main jQuery call
$(document)
    .ready(
	function() {
	    addNavigation();

	    formatTable();

	    initScore();



	}
    );
