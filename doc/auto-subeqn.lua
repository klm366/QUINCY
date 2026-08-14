-- Auto-numbering for subequations using @a, @b, @c placeholders
-- User writes: \tag{@a}, \tag{@b}, \tag{@c}
-- Filter replaces with: \tag{4.29a}, \tag{4.29b}, \tag{4.29c}

local equation_counter = 0
local current_section = 0

-- label -> displayed number ("3.10", "3.10b"), filled in as equations are numbered.
-- References are resolved from this at build time rather than by MathJax, see the
-- resolve_reference pass at the bottom of this file.
local numbers = {}
local anchors_for = {}

-- Cross-file references ------------------------------------------------------
--
-- A filter only ever sees one document, so a reference to an equation on another
-- page cannot be resolved from this render alone. Each render therefore leaves a
-- record of what it numbered in EQREF_DIR, and every render reads the records
-- left by the others. A page split across files resolves after both have been
-- rendered once; a full build renders everything, so it is consistent from then on.
--
-- The directory holds generated data and belongs in .gitignore. A stale record
-- committed to the repository would produce plausible but wrong numbers, which is
-- worse than none at all.
-- Anchored to the project directory, never to the working directory: Quarto renders
-- each file with the working directory set to that file's own folder, so a relative
-- path would scatter a records folder through the whole source tree.
local EQREF_DIR = nil

local source_file = nil   -- this document, as a path below the project directory
local elsewhere = {}      -- label -> { number = "3.10", file = "guide/....qmd" }

local function slug(path)
  return (path:gsub("[^%w]", "_"))
end

local function read_records()
  if not EQREF_DIR then return end
  local listing = io.popen('ls "' .. EQREF_DIR .. '" 2>/dev/null')
  if not listing then return end
  for name in listing:lines() do
    local handle = io.open(EQREF_DIR .. "/" .. name, "r")
    if handle then
      local file = nil
      for line in handle:lines() do
        local path = line:match("^#file\t(.+)$")
        if path then
          file = path
        else
          local label, number = line:match("^(.-)\t(.+)$")
          if label and file then elsewhere[label] = { number = number, file = file } end
        end
      end
      handle:close()
    end
  end
  listing:close()
end

local function write_record()
  if not source_file or not EQREF_DIR then return end
  pandoc.system.make_directory(EQREF_DIR, true)
  local handle = io.open(EQREF_DIR .. "/" .. slug(source_file), "w")
  if not handle then return end
  handle:write("#file\t" .. source_file .. "\n")
  local sorted = {}
  for label in pairs(numbers) do sorted[#sorted + 1] = label end
  table.sort(sorted)
  for _, label in ipairs(sorted) do handle:write(label .. "\t" .. numbers[label] .. "\n") end
  handle:close()
end

-- Drop the group's own \label{} where its rows have been given labels of their own.
-- The rows carry \label{outer_a}, \label{outer_b} and so on, so leaving the group label
-- loose in the same align declares it twice and MathJax reports "Multiple \label".
-- The group stays reachable: record_numbers maps it to the group number and with_anchors
-- emits an anchor for it.
function strip_floating_label(text)
  local lines = {}
  for line in (text .. "\n"):gmatch("([^\n]*)\n") do lines[#lines + 1] = line end
  local kept = {}
  for _, line in ipairs(lines) do
    local lone = line:match("^%s*\\label{([^}]+)}%s*$")
    local drop = false
    if lone then
      local escaped = lone:gsub("([^%w])", "%%%1")
      -- the same label, or a lettered row of it, is attached to a numbered row
      drop = text:match("\\label{" .. escaped .. "_%a}%s*\\tag{") ~= nil
          or text:match("\\label{" .. escaped .. "}%s*\\tag{") ~= nil
    end
    if not drop then kept[#kept + 1] = line end
  end
  return table.concat(kept, "\n")
end

-- Read the label/tag pairs back out of a block once it has been numbered. Working from
-- the finished text keeps this in step with every branch above, rather than repeating
-- their numbering rules.
function record_numbers(text)
  local found = {}
  -- a label immediately followed by its tag: sub-equation rows, and rows carrying a
  -- hand-written label
  for label, tag in text:gmatch("\\label{([^}]+)}%s*\\tag{([^}]+)}") do
    numbers[label] = tag
    found[#found + 1] = label
  end
  -- a plain equation: one label, one tag, anywhere in the block
  local labels, tags = {}, {}
  for l in text:gmatch("\\label{([^}]+)}") do labels[#labels + 1] = l end
  for t in text:gmatch("\\tag{([^}]+)}") do tags[#tags + 1] = t end
  if #labels == 1 and #tags == 1 and not numbers[labels[1]] then
    numbers[labels[1]] = tags[1]
    found[#found + 1] = labels[1]
  end
  -- the group label of a sub-equation block refers to the group as a whole
  local outer = text:match("\\label{([^}]+)_a}")
  local base = text:match("\\tag{(%d+%.%d+)a}")
  if outer and base and not numbers[outer] then
    numbers[outer] = base
    found[#found + 1] = outer
  end
  anchors_for[text] = found
end

-- Emit an anchor of our own for every label in the block. MathJax also creates anchors,
-- but only while typesetting, so with lazy typesetting they do not exist until the
-- reader has scrolled past the equation. These are in the HTML from the start.
function with_anchors(el)
  local labels = anchors_for[el.text]
  if not labels or #labels == 0 then return el end
  local spans = {}
  for _, label in ipairs(labels) do
    spans[#spans + 1] = string.format('<span id="eq:%s" class="eq-anchor"></span>', label)
  end
  return { pandoc.RawInline("html", table.concat(spans)), el }
end

-- Identify this document and pick up the records left by other pages. When the
-- documentation is split across files, each file states which chapter it is with
-- `eq-chapter` in its front matter, and its equations are numbered from there.
-- Without it the chapters are counted within the document, as they always were.
function Meta(meta)
  if quarto and quarto.project and quarto.project.directory
     and quarto.doc and quarto.doc.input_file then
    local root = quarto.project.directory
    EQREF_DIR = root .. "/.eqrefs"
    source_file = quarto.doc.input_file:gsub("^" .. root:gsub("([^%w])", "%%%1") .. "/?", "")
  end

  local chapter = meta["eq-chapter"]
  if chapter then
    local n = tonumber(pandoc.utils.stringify(chapter))
    -- the first level 1 heading takes it to exactly this number
    if n then current_section = n - 1 end
  end

  read_records()
  return meta
end

function Header(el)
  -- Count level 1 headers to track section numbers
  if el.level == 1 then
    current_section = current_section + 1
    equation_counter = 0  -- Reset equation counter for new section
  end
  
  return el
end

function Math(el)
  if el.mathtype ~= "DisplayMath" then
    return el
  end
  
  local text = el.text
  
  -- Check if this is an align environment with an outer label
  local has_align = text:match("\\begin{align}")
  local outer_label = text:match("\\label{([^}]+)}")
  
  -- Auto-add @letter tags to multi-line align blocks with outer labels
  if has_align and outer_label and not text:match("\\tag{") then
    -- Split text into lines
    local lines = {}
    for line in (text .. "\n"):gmatch("([^\n]*)\n") do
      table.insert(lines, line)
    end
    
    -- Find equation lines (those with & and ending with \\)
    local in_align = false
    local letter_index = 0
    local letters = "abcdefghijklmnopqrstuvwxyz"
    local new_lines = {}
    local equation_line_indices = {}
    
    for i, line in ipairs(lines) do
      if line:match("\\begin{align}") then
        in_align = true
        table.insert(new_lines, line)
      elseif line:match("\\end{align}") then
        in_align = false
        table.insert(new_lines, line)
      elseif in_align and line:match("&") then
        -- This is an equation line. Rows that already carry their own \label{}
        -- count too: skipping them used to leave the row untagged, so MathJax's
        -- ams numbering gave it a document-wide number like "(2)" in the middle
        -- of a 6.8a/6.8b/6.8c group, and the remaining letters skipped a step.
        -- The tag replacement below keeps an existing label and only fills in
        -- the letter, so such a row now numbers in sequence with its siblings.
        letter_index = letter_index + 1
        table.insert(equation_line_indices, #new_lines + 1)  -- Remember where to add tag
        table.insert(new_lines, line)
      else
        table.insert(new_lines, line)
      end
    end
    
    -- Add tags based on number of equations found
    if letter_index >= 2 then
      -- Multiple equations - add letter tags
      for idx, line_num in ipairs(equation_line_indices) do
        local line = new_lines[line_num]
        local letter = letters:sub(idx, idx)
        if line:match("\\\\%s*$") then
          new_lines[line_num] = line:gsub("(\\\\%s*)$", "  \\tag{@" .. letter .. "}%1")
        else
          new_lines[line_num] = line .. "  \\tag{@" .. letter .. "}"
        end
      end
      text = table.concat(new_lines, "\n")
    elseif letter_index == 1 then
      -- Single equation - add number tag without letter
      local line_num = equation_line_indices[1]
      local line = new_lines[line_num]
      if line:match("\\\\%s*$") then
        new_lines[line_num] = line:gsub("(\\\\%s*)$", "  \\tag{@0}%1")
      else
        new_lines[line_num] = line .. "  \\tag{@0}"
      end
      text = table.concat(new_lines, "\n")
    end
  end
  
  -- Check if this math block contains subequation tags (@a, @b, @c, etc.)
  if text:match("\\tag{@[a-z]}") or text:match("\\tag{@0}") then
    -- Increment equation counter for this subequation group
    equation_counter = equation_counter + 1
    
    -- Build base number: use section.equation format
    local base_num = string.format("%d.%d", current_section, equation_counter)
    
    -- Extract the outer label if present (e.g., \label{eq:outer})
    local outer_label = text:match("\\label{([^}]+)}")
    
    -- Replace all @letter placeholders with numbered versions
    -- Also add automatic labels if outer_label exists and no label on that line
    text = text:gsub("([^\n]*\\tag{@([a-z0])})([^\n]*)", function(before_tag, letter, after_tag)
      local full_line = before_tag .. after_tag
      local new_tag
      
      if letter == "0" then
        -- Single equation - just use the base number
        new_tag = string.format("\\tag{%s}", base_num)
      else
        -- Multiple equations - add letter suffix
        new_tag = string.format("\\tag{%s%s}", base_num, letter)
      end
      
      -- Check if this line already has a \label{}
      if not full_line:match("\\label{") and outer_label then
        -- No label exists, add one based on outer_label + letter
        local auto_label
        if letter == "0" then
          auto_label = string.format("\\label{%s}", outer_label)
        else
          auto_label = string.format("\\label{%s_%s}", outer_label, letter)
        end
        -- Insert label before the tag
        return full_line:gsub("\\tag{@" .. letter .. "}", auto_label .. new_tag)
      else
        -- Label already exists or no outer label, just replace the tag
        return full_line:gsub("\\tag{@" .. letter .. "}", new_tag)
      end
    end)
    
    el.text = text
  elseif text:match("\\label{eq:") and not text:match("\\tag{") and not text:match("\\begin{align}") then
    -- Regular single equation with label but no tag and no align - add tag and count it
    equation_counter = equation_counter + 1
    local base_num = string.format("%d.%d", current_section, equation_counter)
    
    -- Find the label line and add tag after the equation content
    local lines = {}
    for line in (text .. "\n"):gmatch("([^\n]*)\n") do
      table.insert(lines, line)
    end
    
    -- Find last non-empty line that's not just whitespace
    local last_content_line_idx = nil
    for i = #lines, 1, -1 do
      if lines[i]:match("%S") then
        last_content_line_idx = i
        break
      end
    end
    
    if last_content_line_idx then
      lines[last_content_line_idx] = lines[last_content_line_idx] .. "  \\tag{" .. base_num .. "}"
    end
    
    text = table.concat(lines, "\n")
    el.text = text
  elseif text:match("\\label{eq:") and not text:match("\\tag{") then
    -- Regular equation with label but no tag - just count it (for align blocks already processed)
    equation_counter = equation_counter + 1
  elseif text:match("\\begin{align}") and not text:match("\\tag{") and not text:match("\\label{") then
    -- Unlabeled align block - add numbered tags with letter suffixes for multiple equations
    equation_counter = equation_counter + 1
    local base_num = string.format("%d.%d", current_section, equation_counter)
    
    -- Split text into lines
    local lines = {}
    for line in (text .. "\n"):gmatch("([^\n]*)\n") do
      table.insert(lines, line)
    end
    
    -- Find all equation lines
    local in_align = false
    local equation_line_indices = {}
    local letters = "abcdefghijklmnopqrstuvwxyz"
    
    for i, line in ipairs(lines) do
      if line:match("\\begin{align}") then
        in_align = true
      elseif line:match("\\end{align}") then
        in_align = false
      elseif in_align and (line:match("&") or line:match("=")) and not line:match("\\\\\\\\") then
        -- This is an equation line (has & or =, but not just \\)
        table.insert(equation_line_indices, i)
      end
    end
    
    -- Add tags based on number of equations found
    if #equation_line_indices >= 2 then
      -- Multiple equations - add letter tags
      for idx, line_num in ipairs(equation_line_indices) do
        local line = lines[line_num]
        local letter = letters:sub(idx, idx)
        if line:match("\\\\%s*$") then
          lines[line_num] = line:gsub("(\\\\%s*)$", "  \\tag{" .. base_num .. letter .. "}%1")
        else
          lines[line_num] = line .. "  \\tag{" .. base_num .. letter .. "}"
        end
      end
    elseif #equation_line_indices == 1 then
      -- Single equation - just use base number
      local line_num = equation_line_indices[1]
      local line = lines[line_num]
      if line:match("\\\\%s*$") then
        lines[line_num] = line:gsub("(\\\\%s*)$", "  \\tag{" .. base_num .. "}%1")
      else
        lines[line_num] = line .. "  \\tag{" .. base_num .. "}"
      end
    end
    
    text = table.concat(lines, "\n")
    el.text = text
  elseif not text:match("\\tag{") and not text:match("\\begin{align}") then
    -- Any other display math without tag and without align - add numbered tag
    equation_counter = equation_counter + 1
    local base_num = string.format("%d.%d", current_section, equation_counter)
    
    -- Find last non-empty line and add tag
    local lines = {}
    for line in (text .. "\n"):gmatch("([^\n]*)\n") do
      table.insert(lines, line)
    end
    
    local last_content_line_idx = nil
    for i = #lines, 1, -1 do
      if lines[i]:match("%S") then
        last_content_line_idx = i
        break
      end
    end
    
    if last_content_line_idx then
      lines[last_content_line_idx] = lines[last_content_line_idx] .. "  \\tag{" .. base_num .. "}"
    end
    
    text = table.concat(lines, "\n")
    el.text = text
  end

  el.text = strip_floating_label(el.text)
  record_numbers(el.text)
  return with_anchors(el)
end

-- Third pass: turn each reference into an ordinary link carrying the number as text.
--
-- \Eqref and \ref are MathJax macros resolved in the browser from labels registered while
-- typesetting. That ties a reference to its target having been typeset, which breaks under
-- lazy typesetting and can never reach another page. The numbers are known here, so the
-- reference is written out as plain HTML instead: it needs no MathJax, it is present in
-- the generated file where it can be checked, and it costs nothing to render.
local function resolve_reference(el)
  if el.mathtype ~= "InlineMath" then return nil end
  local kind, label = el.text:match("^%s*\\(%a+){([^}]+)}%s*$")
  if not (kind == "Eqref" or kind == "ref") then return nil end

  local number = numbers[label]
  local target = "#eq:" .. label

  if not number then
    -- Not in this document: try the records left by the other pages. The link is
    -- written to the source file and Quarto rewrites it to the rendered page, so
    -- the depth of either page does not matter here.
    local other = elsewhere[label]
    if other then
      number = other.number
      target = "/" .. other.file .. "#eq:" .. label
    end
  end

  if not number then
    -- Loud on purpose: an unresolved reference used to appear only as "(???)" in the
    -- browser, with nothing in the build output to notice. On a split document this
    -- also fires when the target page has not been rendered yet, which a second
    -- render resolves.
    io.stderr:write("auto-subeqn.lua: unresolved equation reference '" .. label .. "'\n")
    return pandoc.Strong(pandoc.Str("[?" .. label .. "]"))
  end

  local text = (kind == "Eqref") and ("Eq.\u{00A0}" .. number) or number
  return pandoc.Link({pandoc.Str(text)}, target, "", pandoc.Attr("", {"quarto-eq-ref"}))
end

-- Leave this document's numbers where the other pages can find them.
function Pandoc(doc)
  write_record()
  return doc
end

return {
  {
    Meta = Meta,
    Header = Header,
    Math = Math,
    traverse = 'topdown'
  },
  {
    Math = resolve_reference,
    Pandoc = Pandoc
  }
}
