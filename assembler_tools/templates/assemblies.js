toggleAllContigGroups = function (assemblyId) {
    // Determine the state of the "select all" checkbox
    const assemblyElement = document.getElementById(`cg-all-${assemblyId}`)
    // Select all contig checkboxes related to the assembly

    const contigs = assemblyElement.parentElement.parentElement.querySelectorAll('.contig');

    // Set the checked state of each contig checkbox based on the "select all" checkbox
    contigs.forEach(contig => {
        contigName = contig.getAttribute('data-cg')
        toggleContigGroup(contigName, !assemblyElement.checked)
    });
}

toggleContigGroup = function (eventOrElement, isSelected = null) {
    const cg = typeof eventOrElement === 'string' ? eventOrElement : this.getAttribute('data-cg')

    // get all elements that belong to the same contig
    const cgs = Array.from(document.querySelectorAll(`[data-cg="${cg}"]`))

    // if isSelected is not provided, toggle the selection
    if (isSelected === null) {
        isSelected = cgs[0].classList.contains('selected')
    }

    if (isSelected) {
        // select
        cgs.forEach(cg => {
            cg.classList.remove('selected')
            // move selected path elements to the front so that the border is visible
            if (cg.tagName === 'path') cg.parentNode.prepend(cg)
        })
    } else {
        // deselect
        cgs.forEach(cg => {
            cg.classList.add('selected')
            // move deselected path elements to the back so that the border of other paths may be visible
            if (cg.tagName === 'path') cg.parentNode.appendChild(cg)
        })
    }
}

/* Get dendrogram labels */
function getDendrogramInfo() {
    function extractData(elementId) {
        return JSON.parse(document.getElementById(elementId).textContent)
    }

    function extractLabels(axisSelector) {
        const labels = document.querySelectorAll(axisSelector);
        return Array.from(labels).map(label => label.textContent);
    }

    function isIdentical(textArray1, textArray2) {
        if (textArray1.length !== textArray2.length) return false
        for (let i = 0; i < textArray1.length; i++) {
            if (textArray1[i] !== textArray2[i]) return false
        }
        return true;
    }

    const data = extractData('ani-matrix-data')

    const labels1 = extractLabels('#matplotlib\\.axis_1 text')
    const labels2 = extractLabels('#matplotlib\\.axis_2 text')
    if (!isIdentical(labels1, labels2)) {
        console.error('Dendrogram label mismatch detected. The labels on both axes should be identical for accurate representation. Please verify the input data or the label extraction logic.', {
            axis1Labels: labels1,
            axis2Labels: labels2,
            isIdentical: isIdentical(labels1, labels2)
        })
        alert('The labels on the dendrogram are not identical!')
    }

    return [labels1, data]
}

/* Show popover on hover */
function initPopover(dendrogramLabels, dendrogramData) {
    const paths = document.querySelectorAll('#QuadMesh_1 path');
    const popoverMap = new Map();

    getLabels = function (id) {
        const index = typeof id === 'string' ? parseInt(id.split('-')[1]) : id;
        const labelCol = dendrogramLabels[index % dendrogramLabels.length];
        const labelRow = dendrogramLabels[Math.floor(index / dendrogramLabels.length)];
        return [labelCol, labelRow];
    }

    titleFunction = function (element) {
        const labelCol = element.getAttribute('data-label-col')
        const labelRow = element.getAttribute('data-label-row')
        if (labelCol === labelRow) {
            return `${labelCol}`
        } else {
            return `${labelCol} x ${labelRow}`
        }
    }

    contentFunction = function (element) {
        const labelCol = element.getAttribute('data-label-col')
        const labelRow = element.getAttribute('data-label-row')
        const similarity = dendrogramData[labelCol][labelRow]
        if (labelCol === labelRow) {
            return `Similarity: <strong>${similarity}</strong><br>Click to select.`
        } else {
            return `Similarity: <strong>${similarity}</strong>`
        }
    }

    paths.forEach((path, index) => {
        // Set a unique ID for each path to target with the popover
        const [labelCol, labelRow] = getLabels(index)
        path.setAttribute('id', `path-${index}`);
        path.setAttribute('data-label-col', labelCol);
        path.setAttribute('data-label-row', labelRow);

        const popoverInstance = new bootstrap.Popover(path, {
            trigger: 'manual',
            html: true,
            title: titleFunction,
            content: contentFunction,
            container: 'body',
            placement: 'top'
        });

        // Store the popover instance in the map
        popoverMap.set(path, popoverInstance)

        // Show popover on mouseenter
        path.addEventListener('mouseenter', function () {
            const popover = popoverMap.get(this);
            if (popover) popover.show()
        });

        // Hide popover on mouseleave
        path.addEventListener('mouseleave', function () {
            const popover = popoverMap.get(this)
            if (popover) popover.hide()
        })

        // Add click event listener to toggle border
        if (labelCol === labelRow) {
            path.classList.add('contig-group');
            path.setAttribute('data-cg', labelCol);
            path.addEventListener('click', toggleContigGroup);
        }
    })
}

function initGfavizPopover(contigs) {
    document.querySelectorAll('.gfaviz-container').forEach((gfavizContainer) => {
        const assembly = gfavizContainer.getAttribute('data-assembly');
        gfavizContainer.querySelectorAll('svg text').forEach((textElement) => {
            // move them to the front
            textElement.parentNode.parentNode.appendChild(textElement.parentNode);

            // get contig name
            const contigOriginalName = textElement.textContent
            const contigId = `${assembly}@${contigOriginalName}`
            const contigGroupId = contigs[contigId].contig_group

            // set class contig and contig-<contig>
            textElement.classList.add('contig-group');
            textElement.setAttribute('data-cg', contigGroupId)

            // add bootstrap.Popover
            const popoverInstance = new bootstrap.Popover(textElement, {
                trigger: 'manual',
                html: true,
                title: contigGroupId,
                content: 'content',
                container: 'body',
                placement: 'top'
            });

            // Show popover on mouseenter
            textElement.addEventListener('mouseenter', function () {
                const popover = popoverInstance
                if (popover) popover.show()
            });

            // Hide popover on mouseleave
            textElement.addEventListener('mouseleave', function () {
                const popover = popoverInstance
                if (popover) popover.hide()
            })

            textElement.addEventListener('click', toggleContigGroup)
        });
    });
};


document.addEventListener('DOMContentLoaded', function () {
    document.getElementById('ani-clustermap-container')
        .querySelectorAll('#axes_3 > [id^="text_"]')
        .forEach((textElement) => {
            textElement.style.pointerEvents = 'none';
        });
    const [dendrogramLabels, dendrogramData] = getDendrogramInfo();
    initPopover(dendrogramLabels, dendrogramData);
    fetch('assemblies.json').then(response => response.json()).then(assemblies => {
        const contigs = {}
        const contigGroups = {}
        // assemblies is a dictionary of assemblies {'flye': ..., 'canu': ...}
        Object.entries(assemblies).forEach(([assembly, data]) => {
            Object.entries(data.contig_groups).forEach(([contigGroup, data]) => {
                contigGroups[contigGroup] = data
                Object.entries(data.contigs).forEach(([contig, data]) => {
                    contigs[contig] = data
                })
            })
        })
        return [assemblies, contigGroups, contigs]
    }).then(([assemblies, contigGroups, contigs]) => {
        console.log({assemblies: assemblies, contigGroups: contigGroups, contigs: contigs})
        initGfavizPopover(contigs)
    });
});