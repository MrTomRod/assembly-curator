export const loadPaf = function (paf) {
    // Define column names
    paf = "queryID\tqueryLen\tqueryStart\tqueryEnd\tstrand\trefID\trefLen\trefStart\trefEnd\tnumResidueMatches\tlenAln\tmapQ\n" + paf
    // Create a TSV parser
    const tsvParser = d3.dsvFormat("\t")
    // Parse the PAF file with header
    let table = tsvParser.parse(paf).map(d => ({
        ...d,
        queryLen: parseInt(d.queryLen),
        queryStart: parseInt(d.queryStart),
        queryEnd: parseInt(d.queryEnd),
        refLen: parseInt(d.refLen),
        refStart: parseInt(d.refStart),
        refEnd: parseInt(d.refEnd),
        numResidueMatches: parseInt(d.numResidueMatches),
        mapQ: parseInt(d.mapQ),
        lenAln: parseInt(d.lenAln)
    }))

    // if strand is '-', switch queryStart and queryEnd
    table.map(row => {
        if (row.strand === '+') {
            [row.queryStartFlip, row.queryEndFlip] = [row.queryStart, row.queryEnd];
        } else {
            [row.queryStartFlip, row.queryEndFlip] = [row.queryEnd, row.queryStart];
        }
        return row;
    });

    // Find longest alignment per queryID
    const findLongestScaffolds = function () {
        // Get the longest alignment for each query contig (queryID)
        let queryToLongestAlignment = new Map()
        for (const row of table) {
            if (!queryToLongestAlignment.has(row['queryID']) || queryToLongestAlignment.get(row['queryID'])['lenAln'] < row['lenAln']) {
                queryToLongestAlignment.set(row['queryID'], row)
            }
        }
        // sort by refStart of the longest alignment
        return new Map([...queryToLongestAlignment.entries()]
            .sort((a, b) => a[1]['refStart'] - b[1]['refStart']))
    }
    const queryToLongestAlignment = findLongestScaffolds()

    table.map(row => {
        row.isLongestAlignment = queryToLongestAlignment.get(row.queryID).lenAln === row.lenAln
        return row
    })

    // If longest alignment is on the - strand, invert the contig
    table = table.map(row => {
        if (queryToLongestAlignment.get(row.queryID).strand === '-') {
            row.queryStartFlip = row.queryLen - row.queryStartFlip
            row.queryEndFlip = row.queryLen - row.queryEndFlip
            row.isFlipped = true
        } else {
            row.isFlipped = false
        }
        return row;
    });

    const createContigToLen = function (table, id, len) {
        /* This function maps ID to len, and sorts the Map by Len in descending order */
        return new Map(table
            .map(d => [d[id], d[len]])  // map ID -> Len
            .sort((a, b) => b[1] - a[1]) // sort in descending order
        )
    }

    // Function: create Map mapping ID to Offset
    const createContigToOffset = function (table, id, len, contigMap) {
        const contigToOffset = new Map()
        let offset = 0
        for (let [contigID, contigLen] of contigMap) {
            contigToOffset.set(contigID, offset)
            offset += contigLen
        }
        contigToOffset.set('__end__', offset)
        return contigToOffset
    }

    // Ref: sort contigs by length, descending order
    const contigToLenRef = createContigToLen(table, 'refID', 'refLen')
    const contigToOffsetRef = createContigToOffset(table, 'refID', 'refLen', contigToLenRef)

    // Query: sort contigs by refStart of the longest alignment
    const sortByQueryStart = function () {
        // Return a queryID -> queryLen sort by d[len] (descending order)
        const queryToLen = new Map()
        for (const [queryID, d] of queryToLongestAlignment) {
            queryToLen.set(queryID, d['queryLen'])
        }
        return queryToLen
    }


    const contigToLenQry = sortByQueryStart()
    const contigToOffsetQry = createContigToOffset(table, 'queryID', 'queryLen', contigToLenQry)

    // Update table: add columns queryStartOffset, queryEndOffset, refStartOffset, refEndOffset
    table = table.map(d => ({
        ...d,
        queryStartOffset: contigToOffsetQry.get(d.queryID) + d.queryStartFlip,
        queryEndOffset: contigToOffsetQry.get(d.queryID) + d.queryEndFlip,
        refStartOffset: contigToOffsetRef.get(d.refID) + d.refStart,
        refEndOffset: contigToOffsetRef.get(d.refID) + d.refEnd
    }))

    table.contigToLenRef = contigToLenRef
    table.contigToOffsetRef = contigToOffsetRef
    table.contigToLenQry = contigToLenQry
    table.contigToOffsetQry = contigToOffsetQry
    return table
}

export const dotplot = function (
    table,
    dotplotElement,
    {
        markerSize = 3,
        colorFwd = 'blue',
        colorRev = 'green',
        title = 'Dotplot',
        xtitle = 'Reference',
        ytitle = 'Query'
    } = {}) {
    const getColor = function (d) {
        return d.isFlipped === (d.strand === '-') ? colorFwd : colorRev
        // return d.isLongestAlignment ? 'red' : d.strand === '-' ? colorFwd : colorRev
    }

    const colors = table.map(d => getColor(d))
    const marker = {size: markerSize, color: colors}

    // Create a dotplot
    const tracePointsStart = {
        x: table.map(d => d.refStartOffset),
        y: table.map(d => d.queryStartOffset),
        mode: 'markers',
        type: 'scatter',
        marker: marker,
        name: 'Start Points',
        text: table.map(d => `
Query ID: ${d.queryID} (${d.queryStart}:${d.queryEnd})<br>
Ref ID: ${d.refID} (${d.refStart}:${d.refEnd})<br>
Strand: ${d.strand}<br>
Len Aln: ${d.lenAln}`),
        hoverinfo: 'text' // Specify that only the text should be displayed on hover
    };

    const tracePointsEnd = {
        x: table.map(d => d.refEndOffset),
        y: table.map(d => d.queryEndOffset),
        mode: 'markers',
        type: 'scatter',
        marker: marker,
        name: 'End Points',
        text: table.map(d => `
Query ID: ${d.queryID} (${d.queryStart}:${d.queryEnd})<br>
Ref ID: ${d.refID} (${d.refStart}:${d.refEnd})<br>
Strand: ${d.strand}<br>
Len Aln: ${d.lenAln}`),
        hoverinfo: 'text' // Specify that only the text should be displayed on hover
    };

    const traceSegments = table.map(d => ({
        x: [d.refStartOffset, d.refEndOffset],
        y: [d.queryStartOffset, d.queryEndOffset],
        mode: 'lines',
        type: 'scatter',
        line: {color: getColor(d)},
    }));

    // Combine all traces
    const data = [tracePointsStart, tracePointsEnd, ...traceSegments];


    const tickmarksRef = [];
    table.contigToLenRef.forEach((len, key) => {
        // Calculate tick value
        const offset = table.contigToOffsetRef.get(key);
        const tickval = offset + (len / 2);
        // Push ticktext and tickval to the tickmarks array
        tickmarksRef.push({ticktext: key, tickval: tickval});
    });
    const annotations = tickmarksRef.map(mark => ({
        x: mark.tickval,
        y: 0, // Adjust this value as per your chart's y-axis range
        xref: 'x',
        text: mark.ticktext,
        showarrow: false,
        xanchor: 'center',
        yanchor: 'top',
        font: {
            size: 10,
            color: 'grey'
        },
        textangle: -90
    }));

    // calculate horizontal tickmarks
    const tickmarksQry = [];
    table.contigToLenQry.forEach((len, key) => {
        // Calculate tick value
        const offset = table.contigToOffsetQry.get(key);
        const tickval = offset + (len / 2);
        // Push ticktext and tickval to the tickmarks array
        tickmarksQry.push({ticktext: key, tickval: tickval});
    })
    // add to annotations
    annotations.push(...tickmarksQry.map(mark => ({
        x: 0, // Adjust this value as per your chart's x-axis range
        y: mark.tickval,
        xref: 'paper',
        text: mark.ticktext,
        showarrow: false,
        xanchor: 'center',
        yanchor: 'top',
        font: {
            size: 10,
            color: 'grey'
        }
    })))

    const tickvalsRef = Array.from(table.contigToOffsetRef.values())
    const tickvalsQry = Array.from(table.contigToOffsetQry.values())
    const layout = {
        title: title,
        xaxis: {
            title: xtitle,
            tickvals: tickvalsRef,
            ticktext: tickvalsRef.map(() => ''),
        },
        yaxis: {
            title: ytitle,
            tickvals: tickvalsQry,
            ticktext: tickvalsQry.map(() => ''),
        },
        hovermode: 'closest',
        showlegend: false,
        annotations: annotations,
    };

    let lastDimensions = {width: undefined, height: undefined};

    // Assuming updatePlot is a function that updates your plot
    function updatePlot() {
        const currentDimensions = {
            width: dotplotElement.clientWidth,
            height: dotplotElement.clientHeight
        };

        if (currentDimensions.width === lastDimensions.width && currentDimensions.height === lastDimensions.height) return;

        const updatedLayout = {...layout, ...currentDimensions};

        // Use Plotly.react to update the plot with the new layout
        Plotly.react(dotplotElement, data, updatedLayout);

        lastDimensions = currentDimensions;
    }

    updatePlot()

    // Options for the observer (which mutations to observe)
    const config = {attributes: true, childList: false, subtree: false};

    // Create an instance of MutationObserver with the callback
    const observer = new MutationObserver(updatePlot);

    // Start observing the target node for configured mutations
    observer.observe(dotplotElement, config);

}

export const delta2paf = function (delta) {
    delta = delta.split(/\r?\n/);

    let paf = ''

    let rname, qname, rlen, qlen, qs, qe, rs, re, strand, NM, cigar, x, y, seen_gt = false;
    delta.forEach((line, index) => {
        var m
        if ((m = /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/.exec(line)) != null) {
            rname = m[1], qname = m[2], rlen = parseInt(m[3]), qlen = parseInt(m[4]);
            seen_gt = true;
            return;
        }
        if (!seen_gt) return;
        var t = line.split(" ");
        if (t.length == 7) {
            for (var i = 0; i < 5; ++i)
                t[i] = parseInt(t[i]);
            strand = ((t[0] < t[1] && t[2] < t[3]) || (t[0] > t[1] && t[2] > t[3])) ? 1 : -1;
            rs = (t[0] < t[1] ? t[0] : t[1]) - 1;
            re = t[1] > t[0] ? t[1] : t[0];
            qs = (t[2] < t[3] ? t[2] : t[3]) - 1;
            qe = t[3] > t[2] ? t[3] : t[2];
            x = y = 0;
            NM = parseInt(t[4]);
            cigar = [];
        } else if (t.length == 1) {
            var d = parseInt(t[0]);
            if (d == 0) {
                var blen = 0, cigar_str = [];
                if (re - rs - x != qe - qs - y) throw Error("inconsisnt alignment");
                cigar.push((re - rs - x) << 4);
                for (var i = 0; i < cigar.length; ++i) {
                    blen += cigar[i] >> 4;
                    cigar_str.push((cigar[i] >> 4) + "MID".charAt(cigar[i] & 0xf));
                }
                paf += [qname, qlen, qs, qe, strand > 0 ? '+' : '-',
                    rname, rlen, rs, re, blen - NM, blen, 0,
                    "NM:i:" + NM, "cg:Z:" + cigar_str.join("")
                ].join("\t") + "\n"
            } else if (d > 0) {
                var l = d - 1;
                x += l + 1, y += l;
                if (l > 0) cigar.push(l << 4);
                if (cigar.length > 0 && (cigar[cigar.length - 1] & 0xf) == 2)
                    cigar[cigar.length - 1] += 1 << 4;
                else cigar.push(1 << 4 | 2); // deletion
            } else {
                var l = -d - 1;
                x += l, y += l + 1;
                if (l > 0) cigar.push(l << 4);
                if (cigar.length > 0 && (cigar[cigar.length - 1] & 0xf) == 1)
                    cigar[cigar.length - 1] += 1 << 4;
                else cigar.push(1 << 4 | 1); // insertion
            }
        }
    })
    return paf
}


export const assembliesToPafMinimap = async function (ref, qry, urlOrData = 'url', params) {
    if (params === undefined) {
        params = [
            '-t', '1',
            '-k', '28',
            '-N', '1000000',
            '-p', '0.000001',
            '--no-long-join',
        ].join(' ')
    }

    const CLI = await new Aioli(["minimap2/2.22"], {
        printInterleaved: false
    });

    const paths = await CLI.mount([
        {name: "ref.fna", [urlOrData]: ref},
        {name: "qry.fna", [urlOrData]: qry}
    ]);

    // const cmd = `minimap2 -cx asm5 ${paths[0]} ${paths[1]}`
    const cmd = `minimap2 ${params} ref.fna qry.fna`
    console.info('running:', cmd)

    // Extract the output
    const output = await CLI.exec(cmd);

    return [output.stdout, output.stderr, cmd]
}



export const assembliesToPafMummer = async function (ref, qry, urlOrData = 'url', params) {
    if (params === undefined) {
        params = '--mincluster 65'
    }

    const CLI = await new Aioli(["mummer4/nucmer/4.0.0rc1"], {
        printInterleaved: false
    });

    const paths = await CLI.mount([
        {name: "ref.fna", [urlOrData]: ref},
        {name: "qry.fna", [urlOrData]: qry}
    ]);

    // const cmd = `minimap2 -cx asm5 ${paths[0]} ${paths[1]}`
    const cmd = `nucmer ${params} ref.fna qry.fna`
    console.info('running:', cmd)

    // Extract the output
    const output = await CLI.exec(cmd);

    return [await CLI.cat("out.delta"), output.stderr, cmd]
}
