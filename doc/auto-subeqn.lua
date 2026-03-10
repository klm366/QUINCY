-- Auto-numbering for subequations using @a, @b, @c placeholders
-- User writes: \tag{@a}, \tag{@b}, \tag{@c}
-- Filter replaces with: \tag{4.29a}, \tag{4.29b}, \tag{4.29c}

local equation_counter = 0
local current_section = 0

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
      elseif in_align and line:match("&") and not line:match("\\label{") then
        -- This is an equation line
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
  elseif text:match("\\label{eq-") and not text:match("\\tag{") and not text:match("\\begin{align}") then
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
  elseif text:match("\\label{eq-") and not text:match("\\tag{") then
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
  
  return el
end

return {
  {
    Header = Header,
    Math = Math,
    traverse = 'topdown'
  }
}
